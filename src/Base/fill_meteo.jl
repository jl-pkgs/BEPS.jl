function q2RH(q::FT, tem::FT)::FT
  # Vapour pressure in mbar
  ea = 0.46 * q * (tem + 273.16) / 100
  es = 6.1078 * exp((17.269 * tem) / (237.3 + tem))
  clamp(ea / es * 100, 0.0, 100.0)
end

"""
Srad: W m-2
temp: degC
rain: m
wind: m/s
hum: specific humidity, q
"""
function fill_meteo!(meteo::ClimateData,
  rad::FT, tem::FT, pre::FT, wind::FT, hum::FT)

  meteo.Srad = rad
  meteo.temp = tem
  meteo.rain = pre / 1000 # m to mm
  meteo.wind = wind
  meteo.LR = -200.0 #  -200.0 means no measured long-wave radiation, the value will be 
  meteo.rh = q2RH(hum, tem)
end

function fill_meteo!(meteo::ClimateData, d::DataFrame, k::Int=1)
  fill_meteo!(meteo, d.rad[k], d.tem[k], d.pre[k], d.wind[k], d.hum[k])
end
