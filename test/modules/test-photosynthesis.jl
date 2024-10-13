@testset "photosynthesis_jl" begin  
  ea = 0.9879918003910609      # kPa
  gb_w = 0.60                  # m s-1
  VCmax25 = 124.13532190023277 # μmol m-2 s-1

  β_soil = 0.9514687124755321
  g0_w = 0.0175
  g1_w = 8.0

  cii = 266.0
  LH_leaf = 68.04830184404008

  ## 白天
  Tleaf = 16.0
  Tleaf_c = 16.5
  Rs_leaf = 378.15760935295305

  r1 = clang.photosynthesis_c(Tleaf, Rs_leaf, ea, gb_w, VCmax25, β_soil, g0_w, g1_w, cii, Tleaf_c, LH_leaf)
  r2 = photosynthesis_jl(Tleaf, Rs_leaf, ea, gb_w, VCmax25, β_soil, g0_w, g1_w, cii, Tleaf_c, LH_leaf)
  @test all(isapprox.(r1, r2, rtol=1e-7))

  ## 夜间
  Tleaf = 16.0
  Tleaf_c = 16.5
  Rs_leaf = 0.0

  r1 = clang.photosynthesis_c(Tleaf, Rs_leaf, ea, gb_w, VCmax25, β_soil, g0_w, g1_w, cii, Tleaf_c, LH_leaf)
  r2 = photosynthesis_jl(Tleaf, Rs_leaf, ea, gb_w, VCmax25, β_soil, g0_w, g1_w, cii, Tleaf_c, LH_leaf)
  @test all(isapprox.(r1, r2, rtol=1e-8))
end
