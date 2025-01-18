using DocStringExtensions: TYPEDFIELDS
using DataFrames: DataFrame

# include("../SPAC/SPAC.jl")

export s_coszs, lai2, VCmax
export fill_meteo!, meteo_pack_jl
export snow_density

include("Leaf.jl")
include("list.jl")
include("lai2.jl")
include("VCmax.jl")
include("snow_density.jl")
include("helper.jl")
include("BEPS_helper.jl")

"""
    s_coszs(jday::Int, j::Int, lat::Float64, lon::Float64)

# Example
```julia
jday, hour, lat, lon = 20, 12, 20., 120.
s_coszs(jday, hour, lat, lon)
```
"""
function s_coszs(jday::Int, hour::Int, lat::Float64, lon::Float64)
  Delta = 0.006918 - 0.399912 * cos(jday * 2π / 365.0) + 0.070257 * sin(jday * 2π / 365.0) -
          0.006758 * cos(jday * 4π / 365.0) + 0.000907 * sin(jday * 4π / 365.0)
  # delta is the declination angle of sun.

  hr = hour + lon / 15.0  # UTC time
  # hr =j*24.0/RTIMES; # local time
  hr > 24 && (hr = hr - 24)
  hr < 0 && (hr = 24 + hr)

  Lat_arc = π * lat / 180.0
  Hsolar1 = (hr - 12.0) * 2.0 * π / 24.0 # local hour angle in arc.

  # sin(h)
  CosZs = cos(Delta) * cos(Lat_arc) * cos(Hsolar1) + sin(Delta) * sin(Lat_arc)
  return CosZs
end

function meteo_pack_jl(Ta::FT, RH::FT) where {FT<:Real}
  ρₐ::FT = 1.292 # ρ_air, kg/m3
  es::FT = cal_es(Ta)
  ea::FT = es * RH / 100
  VPD::FT = es - ea

  q::FT = ea2q(ea)
  cp::FT = cal_cp(q)

  λ::FT = cal_lambda(Ta)
  Δ::FT = cal_slope(Ta) # slope of es
  γ::FT = 0.066         # kPa/K, 
  # lambda = cal_lambda(Ta) # J kg-1
  # psy = cp * 101.13 / (0.622 * lambda)  
  (; ρₐ, cp, VPD, λ, Δ, γ, es, ea, q)
end
