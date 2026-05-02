using BEPS, DataFrames, Test, Dates

@testset "beps_optimize" begin
  LAI_all = readdlm(path_proj("examples/input/p1_lai.txt"))[:]
  forcing_all = deserialize(path_proj("data/p1_forcing"))

  # 仅取前 48 小时加速测试（CI 每次评估全年耗时太久）
  nhours = 48
  forcing = forcing_all[1:nhours]
  LAI = LAI_all[1:ceil(Int, nhours / 24)]
  dates = DateTime(2010):Hour(1):DateTime(2010, 1, 1) + Hour(nhours - 1)

  VegType, SoilType = 25, 8
  kw = (lon=120.5, lat=30.5, clumping=0.85, fix_snowpack=false)

  # 生成合成观测
  ps_true = ParamBEPS(VegType, SoilType)
  ps_true.r_drainage = 0.55
  ps_true.veg.r_root_decay = 0.90
  Ta = Float64(forcing.Tair[1])
  state_true, _ = setup(ps_true; Ta, Tsoil=2.2, θ0=0.4115, z_snow=0.0)

  _, df_ET_true, _ = beps_modern(forcing, LAI, dates; ps=ps_true, state=state_true, kw...)
  obs_LH = df_ET_true.LH

  # 设置初始猜测
  ps_opt = ParamBEPS(VegType, SoilType)
  ps_opt.r_drainage = 0.30
  ps_opt.veg.r_root_decay = 0.98

  # n_complex=2 将初始集团从 25 降至 10，大幅加速
  model_final, best_rmse = beps_optimize(forcing, LAI, dates, ps_opt, obs_LH;
    col_sim=:LH, kw...,
    maxn=10, kstop=3, f_reltol=0.1, n_complex=2
  )

  @test model_final isa ParamBEPS
  @test best_rmse >= 0
  @test isfinite(best_rmse)
end
