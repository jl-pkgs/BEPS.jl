function evaporation_soil(temp_air, temp_g, rh_air, netRad_g, Gheat_g,
  percent_snow_g::TypeRef, depth_water::TypeRef, depth_snow::TypeRef, mass_water_g::TypeRef, mass_snow_g::TypeRef,
  density_snow, swc_g, porosity_g,
  evapo_soil::TypeRef, evapo_water_g::TypeRef, evapo_snow_g::TypeRef)

  ccall((:evaporation_soil, libbeps), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble, Cdouble,
      Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
    temp_air, temp_g, rh_air, netRad_g, Gheat_g, percent_snow_g, depth_water, depth_snow, mass_water_g, mass_snow_g, density_snow, swc_g, porosity_g, evapo_soil, evapo_water_g, evapo_snow_g)
end

# """
# ```julia
# percent_snow_g, depth_water, depth_snow, mass_water_g, mass_snow_g,
# evapo_soil, evapo_water_g, evapo_snow_g =
#   evaporation_soil(temp_air, temp_g, rh_air, netRad_g, Gheat_g, density_snow, swc_g, porosity_g);
# ```
# """
# function evaporation_soil(temp_air, temp_g, rh_air, netRad_g, Gheat_g,
#   # percent_snow_g::TypeRef, depth_water::TypeRef, depth_snow::TypeRef, mass_water_g::TypeRef, mass_snow_g::TypeRef,
#   mass_water_g::TypeRef, 
#   density_snow, swc_g, porosity_g
#   # evapo_soil::TypeRef, evapo_water_g::TypeRef, evapo_snow_g::TypeRef
# )
#   percent_snow_g = init_dbl()
#   depth_water = init_dbl()
#   depth_snow = init_dbl()
#   # mass_water_g2 = init_dbl()
#   # mass_snow_g = init_dbl()

#   evapo_soil = init_dbl()
#   evapo_water_g = init_dbl()
#   evapo_snow_g = init_dbl()

#   ccall((:evaporation_soil, libbeps), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble, Cdouble,
#       Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
#     temp_air, temp_g, rh_air, netRad_g, Gheat_g,
#     percent_snow_g, depth_water, depth_snow, mass_water_g, mass_snow_g,
#     density_snow, swc_g, porosity_g,
#     evapo_soil, evapo_water_g, evapo_snow_g)

#   percent_snow_g[], depth_water[], depth_snow[],
#   # mass_water_g[], 
#   mass_snow_g[],
#   evapo_soil[], evapo_water_g[], evapo_snow_g[]
# end
