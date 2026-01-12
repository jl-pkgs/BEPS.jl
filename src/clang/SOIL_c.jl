import BEPS: soil_root_fraction!,
  init_soil_parameters!,
  init_soil_status!,
  soil_water_factor_v2,
  Root_Water_Uptake,
  Update_Cs,
  Update_ice_ratio,
  UpdateSoilThermalConductivity,
  UpdateHeatFlux,
  Update_Tsoil_c,
  Update_G,
  UpdateSoilMoisture

function soil_root_fraction!(p::Soil_c)
  ccall((:SoilRootFraction, libbeps), Cvoid, (Ptr{Soil_c},), Ref(p))
end

function init_soil_parameters!(landcover::Int, stxt::Int, r_root_decay::Real, p::Soil_c)
  ccall((:Init_Soil_Parameters, libbeps), Cvoid,
    (Cint, Cint, Cdouble, Ptr{Soil_c}),
    Cint(landcover), Cint(stxt), Cfloat(r_root_decay), Ref(p))
end

function init_soil_status!(p::Soil_c, Tsoil, Tair, Ms, snowdepth)
  ccall((:Init_Soil_Status, libbeps), Cvoid,
    (Ptr{Soil_c}, Cdouble, Cdouble, Cdouble, Cdouble),
    Ref(p), Tsoil, Tair, Ms, snowdepth)
end

function soil_water_factor_v2(p::Soil_c)
  ccall((:soil_water_factor_v2, libbeps), Cvoid, (Ptr{Soil_c},), Ref(p))
end

function Root_Water_Uptake(p::Soil_c, Trans_o, Trans_u, Evap_soil)
  ccall((:Soil_Water_Uptake, libbeps), Cvoid,
    (Ptr{Soil_c}, Cdouble, Cdouble, Cdouble),
    Ref(p), Trans_o, Trans_u, Evap_soil)
end

# function UpdateSoilLambda(p::Soil_c)
#   ccall((:UpdateSoilLambda, libbeps), Cvoid, (Ptr{Soil_c},), Ref(p))
# end

# function init_soil_parameter(T_USDA, S_USDA, Ref_Depth, T_Density, S_Density, T_OC, S_OC, p::Soil_c)
#   ccall((:init_soil_parameter, libbeps), Cvoid,
#     (Cuchar, Cuchar, Cuchar, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Soil_c}),
#     T_USDA, S_USDA, Ref_Depth, T_Density, S_Density, T_OC, S_OC, Ref(p))
# end

function Update_Cs(p::Soil_c)
  ccall((:Update_Cs, libbeps), Cvoid, (Ptr{Soil_c},), Ref(p))
end

function Update_ice_ratio(p::Soil_c)
  ccall((:Update_ice_ratio, libbeps), Cvoid, (Ptr{Soil_c},), Ref(p))
end

function UpdateSoilThermalConductivity(p::Soil_c)
  ccall((:UpdateSoilThermalConductivity, libbeps), Cvoid, (Ptr{Soil_c},), Ref(p))
end

function UpdateHeatFlux(p::Soil_c, Xcs_g, κ_snow, Tsn0,
  Tair_annual_mean, peroid_in_seconds)

  ccall((:UpdateHeatFlux, libbeps), Cvoid,
    (Ptr{Soil_c}, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble),
    Ref(p), Xcs_g, κ_snow, Tsn0, Tair_annual_mean, peroid_in_seconds)
end

function UpdateHeatFlux(p::Soil_c,
  # Xcs_g, κ_snow, Tsn0,
  Tair_annual_mean, peroid_in_seconds)

  Xcs_g, κ_snow, Tsn0 = 0.0, 0.0, 0.0
  ccall((:UpdateHeatFlux, libbeps), Cvoid,
    (Ptr{Soil_c}, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble),
    Ref(p), Xcs_g, κ_snow, Tsn0, Tair_annual_mean, peroid_in_seconds)
end

function Update_Tsoil_c(p::Soil_c, value::Cdouble)
  ccall((:Update_Tsoil_c, libbeps), Cvoid, (Ptr{Soil_c}, Cdouble),
    Ref(p), value)
end

function Update_G(p::Soil_c, value::Cdouble)
  ccall((:Update_G, libbeps), Cvoid, (Ptr{Soil_c}, Cdouble),
    Ref(p), value)
end

function UpdateSoilMoisture(p::Soil_c, peroid_in_seconds::Cdouble)
  ccall((:UpdateSoilMoisture, libbeps), Cvoid, (Ptr{Soil_c}, Cdouble),
    Ref(p), peroid_in_seconds)
end
