# @testset "surface_temperature_jl tests" 
begin
  T_air = 2.0  # Air temperature in Celsius
  rh_air = 0.5  # Relative humidity
  z_snow = 0.03  # Depth of snow in meters
  z_water = 0.1  # Depth of water in meters
  cp_soil1 = 1000.0  # Heat capacity of soil layer 1
  cp_soil0 = 800.0  # Heat capacity of soil layer 0
  Gheat_g = 0.1  # Ground heat flux
  depth_soil1 = 0.5  # Depth of soil layer 1 in meters
  ρ_snow = 300.0  # Density of snow in kg/m^3
  tempL_u = 0.0  # Temperature of the lower layer in Celsius
  Rn_g = 200.0  # Net radiation on the ground
  E_soil = 0.1  # Evaporation from soil
  E_water_g = 0.05  # Evaporation from water on the ground
  E_snow_g = 0.02  # Evaporation from snow on the ground
  λ_soil1 = 0.5  # Thermal conductivity of soil layer 1
  perc_snow_g = 0.3  # Percentage of snow cover on the ground
  G_soil1 = 0.2  # Soil heat flux in layer 1
  T_ground_last = 0.001 # Last ground temperature in Celsius
  T_soil1_last = 0.001 # Last soil layer 1 temperature in Celsius
  T_any0_last = 0.001 # Last temperature of any layer 0 in Celsius
  T_soil0_last = 0.001 # Last soil layer 0 temperature in Celsius
  T_snow_last = 0.001 # Last snow temperature in Celsius
  T_snow1_last = 0.001 # Last snow layer 1 temperature in Celsius
  T_snow2_last = 0.001 # Last snow layer 2 temperature in Celsius

  # T_ground, T_any0, T_soil0, T_snow, T_snow1, T_snow2, G
  for z_snow = [0.01, 0.03, 0.06]
    r_jl = surface_temperature_jl(
      T_air, rh_air,
      z_snow, z_water,
      cp_soil1, cp_soil0, Gheat_g,
      depth_soil1, ρ_snow, tempL_u, Rn_g, E_soil, E_water_g, E_snow_g, λ_soil1, perc_snow_g, G_soil1,
      T_ground_last, T_soil1_last, T_any0_last, T_soil0_last, T_snow_last, T_snow1_last, T_snow2_last
    )
    r_c = clang.surface_temperature_c(
      T_air, rh_air,
      z_snow, z_water,
      cp_soil1, cp_soil0, Gheat_g,
      depth_soil1, ρ_snow, tempL_u, Rn_g, E_soil, E_water_g, E_snow_g, λ_soil1, perc_snow_g, G_soil1,
      T_ground_last, T_soil1_last, T_any0_last, T_soil0_last, T_snow_last, T_snow1_last, T_snow2_last
    )
    # println(r_jl)
    # println(r_c)
    @test all(isapprox.(r_jl, r_c, rtol=1e-8))
  end
end
