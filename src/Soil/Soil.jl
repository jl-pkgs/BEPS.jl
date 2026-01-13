include("SoilHelper.jl")
include("UpdateHeatFlux.jl")
include("UpdateSoilMoisture.jl")
include("soil_water_factor_v2.jl")


function Update_Tsoil_c(p::Soil, value::Cdouble)
  p.Tsoil_c[1] = value
end

function Update_G(p::Soil, value::Cdouble)
  p.G[1] = value
end

export Init_Soil_Parameters, Init_Soil_T_θ!, Sync_Soil_to_State!, UpdateRootFraction, Update_Tsoil_c, Update_G
export UpdateHeatFlux, Update_Cs,
  Update_ice_ratio,
  UpdateSoilThermalConductivity,
  soil_water_factor_v2,
  UpdateSoilMoisture, Root_Water_Uptake
