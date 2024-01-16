using Test
using BEPS

@testset "readparam" begin
  lcs = [1, 2, 6, 9, 13, 40, -1]
  for lc = lcs
    par1 = readparam(lc)
    par2 = clang.readparam(lc)
    @test maximum(abs.(par1 .- par2)) == 0
  end
end

@testset "readcoef" begin
  lcs = [1, 6, 13, -1]
  stxts = 1:11

  for lc = lcs, stxt = stxts
    coef1 = readcoef(lc, stxt)
    coef2 = clang.readcoef(lc, stxt)
    @test maximum(abs.(coef1 .- coef2)) == 0
  end
end

@testset "aerodynamic_conductance" begin
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
