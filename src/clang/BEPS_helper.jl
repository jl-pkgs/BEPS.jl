using DataFrames: DataFrame


function fill_meteo!(meteo::ClimateData, d::DataFrame, k::Int=1)
  meteo.Srad = d.rad[k]
  meteo.temp = d.tem[k]
  meteo.rain = d.pre[k] / 1000 # m to mm
  meteo.wind = d.wind[k]
  meteo.LR = -200.0 #  -200.0 means no measured long-wave radiation, the value will be calculated later

  tem = meteo.temp
  hum = d.hum[k]

  # Vapour pressure in mbar
  es = 0.46 * hum * (tem + 273.16) / 100
  esd = 6.1078 * exp((17.269 * tem) / (237.3 + tem))
  meteo.rh = clamp(es / esd * 100, 0.0, 100.0) # relative humidity, %
  nothing
end


function init_soil!(p_soil::Soil, var_o::Vector,
  meteo::ClimateData, parameter::Vector, par::NamedTuple)
  # landcover::Int, soil_type::Int, Tsoil, soilwater, snowdepth
  # /***** initialize soil conditions, read soil parameters and set depth *****/
  Init_Soil_Parameters(par.landcover, par.soil_type, parameter[28], p_soil)
  p_soil.r_drainage = parameter[27]
  Init_Soil_Status(p_soil, par.Tsoil, meteo.temp, par.soilwater, par.snowdepth) # LHE
  # Initialize intermediate variables array
  var_o .= 0
  for i = 4:9
    var_o[i] = meteo.temp
  end
  for i = 10:15
    var_o[i] = p_soil.temp_soil_p[i-9]
  end
  for i = 22:27
    var_o[i] = p_soil.thetam_prev[i-21]
  end
  for i = 28:33
    var_o[i] = p_soil.ice_ratio[i-27]
  end
  # for (i=0;i<=40;i++)   var_o[i] = 0;
  # for (i=3;i<=8;i++)   var_o[i] = tem;
  # for(i=9;i<=14;i++) var_o[i] = p_soil->temp_soil_p[i-9];
  # for(i=21;i<=26;i++) var_o[i] = p_soil->thetam_prev[i-21];
  # for(i=27;i<=32;i++) var_o[i] = p_soil->ice_ratio[i-27];
  nothing
end


## put struct into a data.frame

function Base.getindex(x::T, i::Int) where {T<:Union{Results,ClimateData}}
  key = fieldnames(T)[i]
  getfield(x, key)
end

Base.length(x::T) where {T<:Union{Results,ClimateData}} = length(fieldnames(T))

function fill_res!(df::DataFrame, Res::T, k::Int) where {T<:Union{Results,ClimateData}}
  n = length(Res)
  for i = 1:n
    df[k, i] = Res[i]
  end
  nothing
end
