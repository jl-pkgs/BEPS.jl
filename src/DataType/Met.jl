export Met

"""
# Fields
$(TYPEDFIELDS)
"""
@with_kw mutable struct Met
  "Inward shortwave radiation, `[W m⁻²]`"
  Srad::Cdouble = 0.0
  
  "Inward longwave radiation, `[W m⁻²]`"
  LR::Cdouble = 0.0

  "2m air temperature, `[°C]`"
  temp::Cdouble = 0.0
  
  "Relative Humidity, `[%]`"
  rh::Cdouble = 0.0
  
  "precipitation, `[mm/h]`"
  rain::Cdouble = 0.0
  
  "2m wind, `[m/s]`"
  wind::Cdouble = 0.0  
end

# Met(Srad, LR, temp, rh, rain, wind) = 
#   Met(; Srad, LR, temp, rh, rain, wind)

"""
# Arguments
- `Srad`: W m-2
- `temp`: degC
- `rain`: mm
- `wind`: m/s
- `hum`: specific humidity, q
"""
function fill_meteo!(met::Met,
  rad::FT, tem::FT, pre::FT, wind::FT, hum::FT)

  met.Srad = rad
  met.temp = tem
  met.rain = pre / 1000 # mm to m
  met.wind = wind
  met.LR = -200.0 #  -200.0 means no measured long-wave radiation, the value will be 
  met.rh = q2RH(hum, tem)
end

function fill_meteo!(met::Met, d::DataFrame, k::Int=1)
  fill_meteo!(met, d.rad[k], d.tem[k], d.pre[k], d.wind[k], d.hum[k])
end


"""
    AirLayer{FT}

Struct to store environmental conditions in each air layer corresponds to one canopy layer.

# Fields
$(TYPEDFIELDS)

> Copied from Land.jl
"""
@with_kw mutable struct AirLayer{FT<:AbstractFloat}
  "Air temperature `[K]`"
  Tair::FT
  "Air density `[kg m⁻³]`"
  ρₐ::FT
  "Specific heat of air `[J kg⁻¹ K⁻¹]`"
  Cp_ca::FT
  "Vapor pressure deficit `[Pa]`"
  VPD::FT
  "Psychrometric constant `[Pa K⁻¹]`"
  γ::FT
  "Slope of Saturation vapor pressure es `[Pa K⁻¹]`"
  Δ::FT
  "Retive humility `[%]`"
  RH::FT
  "Wind speed `[m s⁻¹]`"
  wind::FT = FT(2)
end
