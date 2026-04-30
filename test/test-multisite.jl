using BEPS, Test, DataFrames
using BEPS: path_proj, run_multisite, beps_multisite_optimize, evaluate_multisite

# 复用短时段数据以保持测试速度
const _NHOURS = 48
const _NDAYS  = ceil(Int, _NHOURS / 24)

function _make_data()
  d_all   = deserialize(path_proj("data/p1_meteo"))
  lai_all = readdlm(path_proj("examples/input/p1_lai.txt"))[:]
  d_all[1:_NHOURS, :], lai_all[1:_NDAYS]
end

@testset "run_multisite serial" begin
  d, lai = _make_data()

  sites = DataFrame(
    site_id  = ["sA", "sB"],
    lon      = [120.5, 121.0],
    lat      = [30.5,  31.0],
    VegType  = [25,    25],
    SoilType = [8,     8],
    clumping = [0.85,  0.85],
    Tsoil0   = [2.2,   3.0],
    θ0       = [0.41,  0.38],
    z_snow0  = [0.0,   0.0],
  )
  forcing_dict = Dict("sA" => d, "sB" => d)
  lai_dict     = Dict("sA" => lai, "sB" => lai)

  res = run_multisite(sites, forcing_dict, lai_dict;
          parallel=false, fix_snowpack=false)

  @test res isa Dict
  @test haskey(res, "sA")
  @test haskey(res, "sB")

  for (_, v) in res
    @test v.df_out isa DataFrame
    @test v.df_ET  isa DataFrame
    @test size(v.Tsoil, 1) == _NHOURS
    @test size(v.θ,     1) == _NHOURS
  end
end

@testset "run_multisite missing site warn" begin
  d, lai = _make_data()
  sites = DataFrame(
    site_id=["X"], lon=[120.0], lat=[30.0],
    VegType=[25], SoilType=[8], clumping=[0.85],
    Tsoil0=[2.0], θ0=[0.4], z_snow0=[0.0]
  )
  # only "X" in forcing, missing in lai_dict → @warn is emitted and site is skipped
  res = run_multisite(sites, Dict("X" => d), Dict{String,Vector}();
          parallel=false, fix_snowpack=false)
  @test !haskey(res, "X")   # site skipped due to missing LAI
end

@testset "beps_multisite_optimize basic" begin
  d, lai = _make_data()

  sites = DataFrame(
    site_id  = ["sA", "sB"],
    lon      = [120.5, 120.5],
    lat      = [30.5,  30.5],
    VegType  = [25,    25],
    SoilType = [8,     8],
    clumping = [0.85,  0.85],
    Tsoil0   = [2.2,   2.2],
    θ0       = [0.41,  0.41],
    z_snow0  = [0.0,   0.0],
  )
  forcing_dict = Dict("sA" => d, "sB" => d)
  lai_dict     = Dict("sA" => lai, "sB" => lai)

  # 生成合成观测（使用默认参数运行）
  model_true = ParamBEPS(25, 8)
  df_out_true, _, _, _ = beps_modern(d, lai; model=model_true,
                           lon=120.5, lat=30.5, fix_snowpack=false, verbose=false)
  obs_gpp = df_out_true.GPP
  obs_dict = Dict("sA" => obs_gpp, "sB" => obs_gpp)

  # 从偏离值出发优化
  model_opt = ParamBEPS(25, 8)
  model_opt.veg.VCmax25 *= 1.3

  # 超快速设置：极少迭代，仅验证接口不崩溃
  model_final, best_rmse = beps_multisite_optimize(
    sites, forcing_dict, lai_dict, obs_dict, model_opt;
    col_sim=:GPP,
    opt_paths=[[:veg, :VCmax25]],
    fix_snowpack=false,
    maxn=5, kstop=2, f_reltol=0.1, n_complex=2
  )

  @test model_final isa ParamBEPS
  @test isfinite(best_rmse)
  @test best_rmse >= 0.0
end
