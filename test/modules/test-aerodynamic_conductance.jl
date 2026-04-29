## V2: 独立模块，不依赖 BEPS 导出，方便物理验证
module AC_V2
  const FT = Float64
  include(normpath(joinpath(@__DIR__, "../../src/aerodynamic_conductance_V2.jl")))
end

# ── V1 vs C（仅 Windows 有效） ──────────────────────────────────────
@testset "aerodynamic_conductance V1 vs C" begin
  if !Sys.iswindows()
    @test true skip = true  # libbeps.dll 仅在 Windows 上可用
  else
    canopyh_o = 2.0
    canopyh_u = 0.2
    height_wind_sp = 2.0
    clumping = 0.8
    Ta = 20.0
    wind_sp = 2.0
    GH_o = 100.0
    pai_o = 4.0
    pai_u = 2.0

    ra_o, ra_u, ra_g, Ga_o, Gb_o, Ga_u, Gb_u =
      aerodynamic_conductance_jl(canopyh_o, canopyh_u, height_wind_sp, clumping, Ta, wind_sp, GH_o,
        pai_o, pai_u)
    r1 = (; ra_o, ra_u, ra_g, Ga_o, Gb_o, Ga_u, Gb_u)

    ra_o, ra_u, ra_g, Ga_o, Gb_o, Ga_u, Gb_u =
      clang.aerodynamic_conductance_c(canopyh_o, canopyh_u, height_wind_sp, clumping, Ta, wind_sp, GH_o,
        pai_o, pai_u)
    r2 = (; ra_o, ra_u, ra_g, Ga_o, Gb_o, Ga_u, Gb_u)

    @test maximum(abs.(values(r1) .- values(r2))) <= 1e-8
  end
end

# ── V2 物理正确性测试（跨平台，不需要 C 库）──────────────────────────
@testset "aerodynamic_conductance V2 物理正确性" begin
  # 典型高大森林参数
  h = 20.0;    h_u = 2.0;  z_wind = 40.0
  Ω = 0.8;    Ta = 20.0;  u = 3.0
  lai_o = 4.0; lai_u = 2.0

  # 1. 零风速 → 返回默认值
  ra_o0, ra_u0, ra_g0, Ga_o0, Gb_o0, Ga_u0, Gb_u0 =
    AC_V2.aerodynamic_conductance_jl(h, h_u, z_wind, Ω, Ta, 0.0, 100.0, lai_o, lai_u)
  @test ra_o0 == 0.0
  @test ra_g0 == 300.0
  @test Ga_o0 ≈ 1 / 200.0
  @test Ga_u0 ≈ 1 / 200.0

  # 2. 近中性（H≈0）→ ra_o 接近解析解
  #    u* ≈ u·k/ln((z-d)/z0m), ra_o_neutral ≈ 1/(k·u*)·ln((z-d)/z0h)
  k = 0.4;  d = 0.8h;  z0m = 0.08h;  z0h = 0.1z0m
  ustar_neutral = u * k / log((z_wind - d) / z0m)
  ra_o_expected = 1 / (k * ustar_neutral) * log((z_wind - d) / z0h)

  ra_o_n, ra_u_n, ra_g_n, Ga_o_n, Gb_o_n, Ga_u_n, Gb_u_n =
    AC_V2.aerodynamic_conductance_jl(h, h_u, z_wind, Ω, Ta, u, 1e-6, lai_o, lai_u)
  @test ra_o_n ≈ ra_o_expected rtol = 0.01   # 1% 偏差

  # 3. 物理合理性：稳定 (H<0) > 中性 > 不稳定 (H>0) → ra_o 有序
  ra_o_stable, = AC_V2.aerodynamic_conductance_jl(h, h_u, z_wind, Ω, Ta, u, -200.0, lai_o, lai_u)
  ra_o_unstable, = AC_V2.aerodynamic_conductance_jl(h, h_u, z_wind, Ω, Ta, u, 200.0, lai_o, lai_u)
  @test ra_o_stable >= ra_o_n          # 稳定层结阻力更大
  @test ra_o_unstable <= ra_o_n        # 不稳定层结阻力更小

  # 4. 冠层各层关系：Ga_o > Ga_u（上层冠层更靠近参考高度，阻力更小）
  ra_o, ra_u, ra_g, Ga_o, Gb_o, Ga_u, Gb_u =
    AC_V2.aerodynamic_conductance_jl(h, h_u, z_wind, Ω, Ta, u, 100.0, lai_o, lai_u)
  @test Ga_o > Ga_u                    # 上层传输导度更大
  @test ra_g >= 120                    # 强制最小值
  @test ra_u > 0
  @test Gb_o > 0
  @test Gb_u > 0

  # 5. 叶片边界层导度有限（非零非无穷）
  @test 0 < Gb_o < 1.0
  @test 0 < Gb_u < 1.0
end

# ── V1 vs V2 差异对比（记录已知差异） ───────────────────────────────
@testset "aerodynamic_conductance V1 vs V2 差异对比" begin
  # 使用典型森林参数
  h = 20.0; h_u = 2.0; z_wind = 40.0; Ω = 0.8
  Ta = 20.0; u = 3.0; H = 100.0; lai_o = 4.0; lai_u = 2.0

  ra_o1, ra_u1, ra_g1, Ga_o1, Gb_o1, Ga_u1, Gb_u1 =
    aerodynamic_conductance_jl(h, h_u, z_wind, Ω, Ta, u, H, lai_o, lai_u)
  ra_o2, ra_u2, ra_g2, Ga_o2, Gb_o2, Ga_u2, Gb_u2 =
    AC_V2.aerodynamic_conductance_jl(h, h_u, z_wind, Ω, Ta, u, H, lai_o, lai_u)

  # V2 使用正确 MOST 公式，结果与 V1 存在合理差异
  # 两者均应在合理物理范围内
  for ra in (ra_o1, ra_o2)
    @test 2.0 <= ra <= 500.0
  end
  for ra in (ra_u1, ra_u2)
    @test ra > 0
  end
  @test ra_g1 >= 120
  @test ra_g2 >= 120

  # V2 改进点：使用 z0h 导致 ra_o 偏大（z0h < z0m，对数项更大）
  @test ra_o2 > ra_o1

  # 打印两版本对比（便于手动核查）
  println("\n── V1 vs V2 数值对比 (h=20m, u=3m/s, H=100W/m²) ──")
  println("       │   V1   │   V2")
  println("ra_o   │ $(round(ra_o1; digits=2)) │ $(round(ra_o2; digits=2))")
  println("ra_u   │ $(round(ra_u1; digits=2)) │ $(round(ra_u2; digits=2))")
  println("ra_g   │ $(round(ra_g1; digits=2)) │ $(round(ra_g2; digits=2))")
  println("Ga_o   │ $(round(Ga_o1; digits=4)) │ $(round(Ga_o2; digits=4))")
  println("Gb_o   │ $(round(Gb_o1; digits=4)) │ $(round(Gb_o2; digits=4))")
end
