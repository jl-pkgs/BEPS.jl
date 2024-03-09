export Met

@with_kw mutable struct Met
  Srad::Cdouble = 0.0
  LR::Cdouble = 0.0
  temp::Cdouble = 0.0
  rh::Cdouble = 0.0
  rain::Cdouble = 0.0
  wind::Cdouble = 0.0
  dr_o::Cdouble = 0.0
  df_o::Cdouble = 0.0
  dr_u::Cdouble = 0.0
  df_u::Cdouble = 0.0
end

Met(Srad, LR, temp, rh, rain, wind) = 
  Met(; Srad, LR, temp, rh, rain, wind)


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
