using Test, BEPS

@testset "netRadiation_jl" begin
  Rs_global = 500.0  # Global solar radiation
  CosZs = 0.5  # Cosine of the solar zenith angle
  temp_o = 300.0  # Overstory temperature in Kelvin
  temp_u = 295.0  # Understory temperature in Kelvin
  temp_g = 290.0  # Ground temperature in Kelvin
  lai_o = 2.0  # Leaf area index of overstory
  lai_u = 1.0  # Leaf area index of understory
  lai_os = 1.5  # Leaf area index of overstory snow
  lai_us = 0.5  # Leaf area index of understory snow
  lai = Leaf(2.0)  # Leaf area index
  Ω = 0.5  # Clumping index
  Tair = 298.0  # Air temperature in Kelvin
  RH = 0.6  # Relative humidity
  α_snow_v = 0.8  # Albedo of snow in visible spectrum
  α_snow_n = 0.7  # Albedo of snow in near-infrared spectrum
  perc_snow_o = 0.3  # Percentage area of snow in overstory
  perc_snow_u = 0.2  # Percentage area of snow in understory
  perc_snow_g = 0.1  # Percentage area of snow on ground
  α_v = Layer3(0.2, 0.25, 0.15)
  α_n = Layer3(0.3, 0.35, 0.25)
  
  Rn_Leaf = Leaf(0.0)  # Net radiation for leaf
  Rns_Leaf = Leaf(0.0)  # Net shortwave radiation for leaf
  Rnl_Leaf = Leaf(0.0)  # Net longwave radiation for leaf
  Ra = Radiation()  # Radiation object

  # Call the function with the test inputs
  r_jl = netRadiation_jl(Rs_global, CosZs, 
    temp_o, temp_u, temp_g, 
    lai_o, lai_u, lai_os, lai_us, lai, Ω, 
    Tair, RH, 
    α_snow_v, α_snow_n, α_v, α_n,
    perc_snow_o, perc_snow_u, perc_snow_g, 
    
    Rn_Leaf, Rns_Leaf, Rnl_Leaf, Ra)
  r_c = clang.netRadiation_c(Rs_global, CosZs,
    temp_o, temp_u, temp_g,
    lai_o, lai_u, lai_os, lai_us, lai, Ω,
    Tair, RH,
    α_snow_v, α_snow_n, α_v, α_n,
    perc_snow_o, perc_snow_u, perc_snow_g,
    Rn_Leaf, Rns_Leaf, Rnl_Leaf, Ra)
  @test all(isapprox.(r_jl, r_c, rtol=1e-8))
end
