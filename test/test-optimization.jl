using BEPS, DataFrames, Test
using BEPS: path_proj

@testset "beps_optimize" begin
  lai = readdlm(path_proj("examples/input/p1_lai.txt"))[:]
  d = deserialize(path_proj("data/p1_meteo"))
  d.tem = d.tem .- 5.0

  kw = (lon=120.5, lat=30.5,
    VegType=25, SoilType=8,
    clumping=0.85,
    Tsoil0=2.2, θ0=0.4115, z_snow0=0.0
  )

  # 用"真值"模型生成合成观测
  model_true = ParamBEPS(kw.VegType, kw.SoilType)
  model_true.r_drainage = 0.55
  model_true.veg.r_root_decay = 0.90

  _, df_ET_true, _, _ = beps_modern(d, lai; model=model_true, kw...,
    fix_snowpack=false, verbose=false)
  obs_LH = df_ET_true.LH

  # 从偏离的初始猜测开始优化
  model_opt = ParamBEPS(kw.VegType, kw.SoilType)
  model_opt.r_drainage = 0.30
  model_opt.veg.r_root_decay = 0.98

  model_final, best_rmse = beps_optimize(d, lai, model_opt, obs_LH;
    col_sim=:LH, kw...,
    fix_snowpack=false, verbose=false,
    maxn=10, kstop=3, f_reltol=0.1
  )

  @test model_final isa ParamBEPS
  @test best_rmse >= 0
  @test isfinite(best_rmse)
end
