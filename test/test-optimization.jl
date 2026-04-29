using BEPS, DataFrames, Test
using BEPS: path_proj

@testset "beps_optimize" begin
  lai_all = readdlm(path_proj("examples/input/p1_lai.txt"))[:]
  d_all = deserialize(path_proj("data/p1_meteo"))
  d_all.tem = d_all.tem .- 5.0

  # 仅取前 48 小时加速测试（CI 每次评估全年耗时太久）
  nhours = 48
  d = d_all[1:nhours, :]
  lai = lai_all[1:ceil(Int, nhours / 24)]

  kw = (lon=120.5, lat=30.5,
    VegType=25, SoilType=8,
    clumping=0.85,
    Tsoil0=2.2, θ0=0.4115, z_snow0=0.0
  )

  # 生成合成观测
  model_true = ParamBEPS(kw.VegType, kw.SoilType)
  model_true.r_drainage = 0.55
  model_true.veg.r_root_decay = 0.90

  _, df_ET_true, _, _ = beps_modern(d, lai; model=model_true, kw...,
    fix_snowpack=false, verbose=false)
  obs_LH = df_ET_true.LH

  # 设置初始猜测
  model_opt = ParamBEPS(kw.VegType, kw.SoilType)
  model_opt.r_drainage = 0.30
  model_opt.veg.r_root_decay = 0.98

  # n_complex=2 将初始集团从 25 降至 10，大幅加速
  model_final, best_rmse = beps_optimize(d, lai, model_opt, obs_LH;
    col_sim=:LH, kw...,
    fix_snowpack=false, verbose=false,
    maxn=10, kstop=3, f_reltol=0.1, n_complex=2
  )

  @test model_final isa ParamBEPS
  @test best_rmse >= 0
  @test isfinite(best_rmse)
end

@testset "beps_optimize with opt_paths" begin
  lai_all = readdlm(path_proj("examples/input/p1_lai.txt"))[:]
  d_all = deserialize(path_proj("data/p1_meteo"))

  nhours = 48
  d   = d_all[1:nhours, :]
  lai = lai_all[1:ceil(Int, nhours / 24)]
  kw  = (lon=120.5, lat=30.5, VegType=25, SoilType=8, clumping=0.85,
         Tsoil0=2.2, θ0=0.4115, z_snow0=0.0)

  model_true = ParamBEPS(kw.VegType, kw.SoilType)
  df_out_ref, _, _, _ = beps_modern(d, lai; model=model_true, kw...,
                          fix_snowpack=false, verbose=false)
  obs_GPP = df_out_ref.GPP

  model_opt = ParamBEPS(kw.VegType, kw.SoilType)
  model_opt.veg.VCmax25 *= 1.5  # 偏离真实值

  # 只优化 VCmax25，其他参数固定
  model_final, best_rmse = beps_optimize(d, lai, model_opt, obs_GPP;
    col_sim=:GPP,
    opt_paths=[[:veg, :VCmax25]],
    kw..., fix_snowpack=false, verbose=false,
    maxn=10, kstop=3, f_reltol=0.1, n_complex=2
  )

  @test model_final isa ParamBEPS
  @test isfinite(best_rmse)
end

@testset "beps_optimize multi-objective" begin
  lai_all = readdlm(path_proj("examples/input/p1_lai.txt"))[:]
  d_all = deserialize(path_proj("data/p1_meteo"))

  nhours = 48
  d   = d_all[1:nhours, :]
  lai = lai_all[1:ceil(Int, nhours / 24)]
  kw  = (lon=120.5, lat=30.5, VegType=25, SoilType=8, clumping=0.85,
         Tsoil0=2.2, θ0=0.4115, z_snow0=0.0)

  model_true = ParamBEPS(kw.VegType, kw.SoilType)
  df_out_ref, df_ET_ref, _, _ = beps_modern(d, lai; model=model_true, kw...,
                                  fix_snowpack=false, verbose=false)

  # 多目标：GPP + LH
  obs = Dict(:GPP => df_out_ref.GPP, :LH => df_ET_ref.LH)

  model_opt = ParamBEPS(kw.VegType, kw.SoilType)
  model_opt.veg.VCmax25 *= 1.3

  model_final, best_f = beps_optimize(d, lai, model_opt, obs;
    kw..., fix_snowpack=false, verbose=false,
    maxn=5, kstop=2, f_reltol=0.1, n_complex=2
  )

  @test model_final isa ParamBEPS
  @test isfinite(best_f)
  @test best_f >= 0.0
end
