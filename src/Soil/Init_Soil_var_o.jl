function Init_Soil_var_o(p_soil::AbstractSoil, var_o::Vector,
  meteo::Met, parameter::Vector, par::NamedTuple)
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
    var_o[i] = p_soil.Tsoil_p[i-9]
  end
  for i = 22:27
    var_o[i] = p_soil.θ_prev[i-21]
  end
  for i = 28:33
    var_o[i] = p_soil.ice_ratio[i-27]
  end
  # for (i=0;i<=40;i++)   var_o[i] = 0;
  # for (i=3;i<=8;i++)   var_o[i] = tem;
  # for(i=9;i<=14;i++) var_o[i] = p_soil->Tsoil_p[i-9];
  # for(i=21;i<=26;i++) var_o[i] = p_soil->θ_prev[i-21];
  # for(i=27;i<=32;i++) var_o[i] = p_soil->ice_ratio[i-27];
  nothing
end
