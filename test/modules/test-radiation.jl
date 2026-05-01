using Test, BEPS

function leaf_values(x::Leaf)
  (x.o_sunlit, x.o_shaded, x.u_sunlit, x.u_shaded)
end

function test_net_radiation_case(; Rs_global, CosZs, T, lai_o, lai_u,
  lai_os, lai_us, lai, Ω, Tair, RH, α_snow_v, α_snow_n,
  perc_snow_o, perc_snow_u, perc_snow_g, α_v, α_n)

  Rn_jl = Leaf()
  Rns_jl = Leaf()
  Rnl_jl = Leaf()
  Ra_jl = Radiation()
  Rn_c = Leaf()
  Rns_c = Leaf()
  Rnl_c = Leaf()
  Ra_c = Radiation()

  r_jl = netRadiation_jl(Rs_global, CosZs, T,
    lai_o, lai_u, lai_os, lai_us, lai, Ω,
    Tair, RH, NaN,
    α_snow_v, α_snow_n, α_v, α_n,
    perc_snow_o, perc_snow_u, perc_snow_g,
    Rn_jl, Rns_jl, Rnl_jl, Ra_jl)
  r_c = clang.netRadiation_c(Rs_global, CosZs, T,
    lai_o, lai_u, lai_os, lai_us, lai, Ω,
    Tair, RH,
    α_snow_v, α_snow_n, α_v, α_n,
    perc_snow_o, perc_snow_u, perc_snow_g,
    Rn_c, Rns_c, Rnl_c, Ra_c)

  @test all(isapprox.(r_jl, r_c; rtol=1e-12, atol=1e-10))
  @test all(isapprox.(leaf_values(Rn_jl), leaf_values(Rn_c); rtol=1e-12, atol=1e-10))
  @test all(isapprox.(leaf_values(Rns_jl), leaf_values(Rns_c); rtol=1e-12, atol=1e-10))
end

@testset "netRadiation_jl" begin
  test_net_radiation_case(
    Rs_global = 500.0,
    CosZs = 0.5,
    T = Layer3(23.0, 22.0, 21.5),
    lai_o = 2.0,
    lai_u = 1.0,
    lai_os = 2.2,
    lai_us = 1.1,
    lai = Leaf(1.1, 0.9, 0.7, 0.4),
    Ω = 0.5,
    Tair = 23.6,
    RH = 65.0,
    α_snow_v = 0.8,
    α_snow_n = 0.7,
    perc_snow_o = 0.3,
    perc_snow_u = 0.2,
    perc_snow_g = 0.1,
    α_v = Layer3(0.2, 0.25, 0.15),
    α_n = Layer3(0.3, 0.35, 0.25),
  )

  test_net_radiation_case(
    Rs_global = 385.5,
    CosZs = 0.35,
    T = Layer3(23.6, 23.6, 23.6),
    lai_o = 2.16,
    lai_u = 0.48,
    lai_os = 2.28,
    lai_us = 0.58,
    lai = Leaf(0.7, 1.58, 0.08, 0.5),
    Ω = 0.85,
    Tair = 23.6,
    RH = 65.61465073124307,
    α_snow_v = 0.0,
    α_snow_n = 0.0,
    perc_snow_o = 0.0,
    perc_snow_u = 0.0,
    perc_snow_g = 0.0,
    α_v = Layer3(0.04, 0.04, 0.1),
    α_n = Layer3(0.42, 0.42, 0.2),
  )
end
