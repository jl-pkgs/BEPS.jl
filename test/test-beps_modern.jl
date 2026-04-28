using BEPS, DataFrames, Test
using BEPS: path_proj

lai = readdlm(path_proj("examples/input/p1_lai.txt"))[:]
kw = (lon=120.5, lat=30.5,
  VegType=25, SoilType=8,
  clumping=0.85,
  Tsoil0=2.2, θ0=0.4115, z_snow0=0.0
)

@testset "beps_modern" begin
  d = deserialize(path_proj("data/p1_meteo"))
  d.tem = d.tem .- 5.0

  # --- 1. 默认参数与 besp_main 结果一致 ---
  @testset "match besp_main (model=nothing)" begin
    df_main, df_ET_main, _, _ = besp_main(d, lai; kw..., version="julia",
      fix_snowpack=false, verbose=false)
    df_mod, df_ET_mod, Tsoil_mod, θ_mod = beps_modern(d, lai; kw...,
      fix_snowpack=false, verbose=false)

    # 输出形状一致
    @test size(df_mod) == size(df_main)
    @test size(df_ET_mod) == size(df_ET_main)

    # 主要通量吻合（容差 1e-7）
    @test maximum(abs.(Matrix(df_mod) .- Matrix(df_main))) <= 1e-7
    @test maximum(abs.(Matrix(df_ET_mod) .- Matrix(df_ET_main))) <= 1e-7

    # 土壤状态的物理边界
    @test all(isfinite, Tsoil_mod)
    @test all(isfinite, θ_mod)
    @test all(θ_mod .>= 0)
    @test all(θ_mod .<= 1)
  end

  # --- 2. 传入自定义 model 时，参数变化能影响输出 ---
  @testset "custom model propagates parameters" begin
    # 基准运行
    df_ref, _, _, _ = beps_modern(d, lai; kw..., fix_snowpack=false, verbose=false)

    # 改变 VCmax25 后输出应有差异
    ps_mod = ParamBEPS(kw.VegType, kw.SoilType)
    ps_mod.veg.VCmax25 *= 2.0
    df_vcmax, _, _, _ = beps_modern(d, lai; model=ps_mod, kw...,
      fix_snowpack=false, verbose=false)

    @test size(df_vcmax) == size(df_ref)
    # GPP 应随 VCmax25 增大而增大
    @test sum(df_vcmax.GPP) > sum(df_ref.GPP)
  end

  # --- 3. 传入 model 时的输出形状 ---
  @testset "explicit model returns correct shape" begin
    ps = ParamBEPS(kw.VegType, kw.SoilType)
    df_m, df_ET_m, Tsoil_m, θ_m = beps_modern(d, lai; model=ps, kw...,
      fix_snowpack=false, verbose=false)

    n = size(d, 1)
    @test size(df_m, 1) == n
    @test size(df_ET_m, 1) == n
    @test size(Tsoil_m) == (n, 5)
    @test size(θ_m) == (n, 5)
  end
end
