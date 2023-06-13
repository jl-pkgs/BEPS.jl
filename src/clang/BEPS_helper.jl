using DataFrames: DataFrame


function hum2RH(hum::FT, tem::FT)::FT
  # Vapour pressure in mbar
  ea = 0.46 * hum * (tem + 273.16) / 100
  es = 6.1078 * exp((17.269 * tem) / (237.3 + tem))
  clamp(ea / es * 100, 0.0, 100.0)
end

function fill_meteo!(meteo::ClimateData,
  rad::FT, tem::FT, pre::FT, wind::FT, hum::FT)

  meteo.Srad = rad
  meteo.temp = tem
  meteo.rain = pre / 1000 # m to mm
  meteo.wind = wind
  meteo.LR = -200.0 #  -200.0 means no measured long-wave radiation, the value will be 
  meteo.rh = hum2RH(hum, tem)
end

function fill_meteo!(meteo::ClimateData, d::DataFrame, k::Int=1)
  fill_meteo!(meteo, d.rad[k], d.tem[k], d.pre[k], d.wind[k], d.hum[k])
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


const TypeDF = Union{Results,ClimateData,OutputET}

## put struct into a data.frame
function Base.getindex(x::T, i::Int)::FT where {T<:TypeDF}
  # key = fieldnames(T)[i]
  getfield(x, i)
end

Base.length(x::T) where {T<:TypeDF} = fieldcount(T)

function fill_res!(df::DataFrame, Res::T, k::Int) where {T<:TypeDF}
  n = length(Res)
  for i in 1:n
    df[k, i] = Res[i]
  end
  nothing
end
