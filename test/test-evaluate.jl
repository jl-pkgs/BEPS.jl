using BEPS, Test, DataFrames
using BEPS: of_NSE, of_KGE, of_R2, of_Bias, of_MAE, evaluate_site, evaluate_multisite

@testset "evaluate metrics" begin
  # 完美预测
  obs = [1.0, 2.0, 3.0, 4.0, 5.0]
  sim = copy(obs)

  @test of_NSE(obs, sim)  ≈ 1.0
  @test of_KGE(obs, sim)  ≈ 1.0
  @test of_R2(obs, sim)   ≈ 1.0
  @test of_Bias(obs, sim) ≈ 0.0 atol=1e-10
  @test of_MAE(obs, sim)  ≈ 0.0 atol=1e-10
  @test of_RMSE(obs, sim) ≈ 0.0 atol=1e-10

  # 系统偏高
  sim2 = obs .+ 1.0
  @test of_NSE(obs, sim2)  < 1.0
  @test of_Bias(obs, sim2) > 0.0   # positive = overestimate
  @test of_MAE(obs, sim2)  ≈ 1.0

  # 完全负相关：R² = cor² = 1（符号被平方掉了）
  obs3 = [1.0, 2.0, 3.0]
  sim3 = [3.0, 2.0, 1.0]  # perfect negative correlation
  @test of_R2(obs3, sim3) ≈ 1.0   # R² = cor²，与相关方向无关

  # NSE < 1 when sim != obs
  @test of_NSE(obs, sim2) < 1.0

  # NaN handling (of_NSE/of_KGE/of_R2/of_Bias/of_MAE handle NaN, of_RMSE is from ModelParams)
  obs_nan = [1.0, NaN, 3.0, 4.0]
  sim_nan = [1.0, 2.0, NaN, 4.0]
  # only indices 1 and 4 valid → perfect match
  @test of_NSE(obs_nan, sim_nan)  ≈ 1.0
  @test of_MAE(obs_nan, sim_nan)  ≈ 0.0 atol=1e-10
end

@testset "evaluate_site returns correct keys" begin
  obs = randn(100) .+ 5
  sim = obs .+ 0.1 .* randn(100)
  m   = evaluate_site(obs, sim)

  @test haskey(m, :RMSE)
  @test haskey(m, :NSE)
  @test haskey(m, :KGE)
  @test haskey(m, :R2)
  @test haskey(m, :Bias)
  @test haskey(m, :MAE)
  @test all(isfinite, values(m))
  @test m.NSE ≤ 1.0
  @test m.KGE ≤ 1.0
  @test 0.0 ≤ m.R2 ≤ 1.0
end

@testset "evaluate_multisite returns DataFrame" begin
  using BEPS: path_proj

  d   = deserialize(path_proj("data/p1_meteo"))
  lai = readdlm(path_proj("examples/input/p1_lai.txt"))[:]
  kw  = (lon=120.5, lat=30.5, VegType=25, SoilType=8, clumping=0.85,
         Tsoil0=2.2, θ0=0.4115, z_snow0=0.0, fix_snowpack=false, verbose=false)

  nhours = 48
  d_sub  = d[1:nhours, :]
  lai_sub = lai[1:ceil(Int, nhours / 24)]

  df_out, df_ET, _, _ = beps_modern(d_sub, lai_sub; kw...)

  # 用模拟值本身作为"观测"（完美精度基准）
  results = Dict("site1" => (; df_out, df_ET, Tsoil=zeros(nhours, 5), θ=zeros(nhours, 5)))
  obs_gpp = df_out.GPP

  df_eval = evaluate_multisite(results, Dict("site1" => obs_gpp); col=:GPP)

  @test df_eval isa DataFrame
  @test nrow(df_eval) == 1
  @test df_eval.site[1] == "site1"
  @test df_eval.RMSE[1] ≈ 0.0 atol=1e-10
  @test df_eval.NSE[1]  ≈ 1.0 atol=1e-8

  # missing site_id should be silently skipped
  df_eval2 = evaluate_multisite(results, Dict("no_such_site" => obs_gpp); col=:GPP)
  @test nrow(df_eval2) == 0
end
