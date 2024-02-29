include("Init_Soil_Parameters.jl")
include("Update_thermal.jl")
include("UpdateSoilMoisture.jl")


function Update_Tsoil_c(p::Soil, value::Cdouble)
  p.Tsoil_c[1] = value
end

function Update_G(p::Soil, value::Cdouble)
  p.G[1] = value
end


export Init_Soil_Parameters, Init_Soil_Status, SoilRootFraction, Update_Tsoil_c, Update_G
export UpdateHeatFlux, Update_Cs,
  Update_ice_ratio,
  UpdateSoilThermalConductivity,
  soil_water_factor_v2,
  UpdateSoilMoisture,
  Soil_Water_Uptake
