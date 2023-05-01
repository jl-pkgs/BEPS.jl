function surface_temperature(temp_air, rh_air, depth_snow, depth_water, capacity_heat_soil1, capacity_heat_soil0, Gheat_g, depth_soil1, density_snow, tempL_u, netRad_g, evapo_soil, evapo_water_g, evapo_snow_g, lambda_soil1, percent_snow_g, heat_flux_soil1, 
  temp_ground_last, temp_soil1_last, temp_any0_last, temp_snow_last, temp_soil0_last, temp_snow1_last, temp_snow2_last,
  temp_ground::TypeRef, temp_any0::TypeRef, temp_snow::TypeRef, temp_soil0::TypeRef, 
  temp_snow1::TypeRef, temp_snow2::TypeRef, heat_flux::TypeRef)
  
  ccall((:surface_temperature, libbeps), Cvoid,
    (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble,
      Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
    temp_air, rh_air, depth_snow, depth_water, capacity_heat_soil1, capacity_heat_soil0, Gheat_g, depth_soil1, density_snow, tempL_u, netRad_g, evapo_soil, evapo_water_g, evapo_snow_g, lambda_soil1, percent_snow_g, heat_flux_soil1, temp_ground_last, temp_soil1_last, temp_any0_last, temp_snow_last, temp_soil0_last, temp_snow1_last, temp_snow2_last,
    temp_ground, temp_any0, temp_snow, temp_soil0, temp_snow1, temp_snow2, heat_flux)
end


function surface_temperature(temp_air, rh_air, depth_snow, depth_water, capacity_heat_soil1, capacity_heat_soil0, Gheat_g, depth_soil1, density_snow, tempL_u, netRad_g, evapo_soil, evapo_water_g, evapo_snow_g, lambda_soil1, percent_snow_g, heat_flux_soil1, temp_ground_last, temp_soil1_last, temp_any0_last, temp_snow_last, temp_soil0_last, temp_snow1_last, temp_snow2_last)

  temp_ground = Ref(0.0)
  temp_any0 = Ref(0.0)
  temp_snow = Ref(0.0)
  temp_soil0 = Ref(0.0)
  temp_snow1 = Ref(0.0)
  temp_snow2 = Ref(0.0)
  heat_flux = Ref(0.0)

  ccall((:surface_temperature, libbeps), Cvoid,
    (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble,
      Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
    temp_air, rh_air, depth_snow, depth_water, capacity_heat_soil1, capacity_heat_soil0, Gheat_g, depth_soil1, density_snow, tempL_u, netRad_g, evapo_soil, evapo_water_g, evapo_snow_g, lambda_soil1, percent_snow_g, heat_flux_soil1, temp_ground_last, temp_soil1_last, temp_any0_last, temp_snow_last, temp_soil0_last, temp_snow1_last, temp_snow2_last, temp_ground, temp_any0, temp_snow, temp_soil0, temp_snow1, temp_snow2, heat_flux)

  temp_ground[], temp_any0[], temp_snow[], temp_soil0[], temp_snow1[], temp_snow2[], heat_flux[]
end
