@testset "Layer3 and Layer2" begin
  x = Layer3()
  y = Layer3(1.0, 2.0)
  @test y.g == 0.0
  z = Layer3(y)

  set!(z, 1.0)
  @test z.g == 1.0

  x = Layer2()
  y = Layer2(1.0)
  @test y.u == 1.0

  z = Layer2(y)
  set!(z, 1.0)
  @test z.o == 1.0
  @test z.u == 1.0

  set!(x, y)
  @test x.o == 1.0
  @test x.u == 1.0
  print(x)
end


@testset "VCmax" begin
  VCmax25 = 62.5          # maximum capacity of Rubisco at 25C-Vcmax	
  N_leaf = 3.10 + 1.35    # leaf Nitrogen content	mean value + 1 SD g/m2 
  slope = 20.72 / 62.5    # slope of Vcmax-N curve
  lai = 1.0
  Ω = 0.5
  CosZs = 0.5
  VCmax_sunlit, VCmax_shaded = VCmax(lai, Ω, CosZs, VCmax25, N_leaf, slope)

  @test VCmax_sunlit == 80.65125251476604
  @test VCmax_shaded == 75.99251585043879
end


@testset "snow_density" begin
  Ta = -50.0:50
  ρ_snow = snow_density.(Ta, 2.0)
  @test minimum(ρ_snow) == 137.33588449248248
  @test maximum(ρ_snow) == 256.4936370728328
end
