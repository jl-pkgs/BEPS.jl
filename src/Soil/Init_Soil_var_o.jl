# landcover::Int, soil_type::Int, Tsoil, soilwater, snowdepth
# initialize soil conditions, read soil parameters and set depth
"Init_Soil_var_o(soil, state, Ta, parameter, par)"
function Init_Soil_var_o(soil::AbstractSoil, state::Union{State,Vector},
  Ta::Real, parameter::Vector, par::NamedTuple)

  Init_Soil_Parameters(par.landcover, par.soil_type, parameter[28], soil)
  soil.r_drainage = parameter[27]
  Init_Soil_Status(soil, par.Tsoil, Ta, par.soilwater, par.snowdepth) # LHE

  init_state!(state, soil, Ta)
end

# for (i=3;i<=8;i++)   var_o[i] = tem;
# for(i=9;i<=14;i++) var_o[i] = soil->Tsoil_p[i-9];
# for(i=21;i<=26;i++) var_o[i] = soil->θ_prev[i-21];
# for(i=27;i<=32;i++) var_o[i] = soil->ice_ratio[i-27];
function init_state!(state::Vector, soil::AbstractSoil, Ta)
  state .= 0
  for i = 1:6
    state[i+3] = Ta
    state[i+9] = soil.Tsoil_p[i]
    state[i+21] = soil.θ_prev[i]
    state[i+27] = soil.ice_ratio[i]
  end
  return nothing
end

function init_state!(state::State, soil::AbstractSoil, Ta)
  state.Ts .= Ta
  state.Tsoil_prev .= soil.Tsoil_p[1:6]
  state.θ_prev .= soil.θ_prev[1:6]
  state.ice_ratio .= soil.ice_ratio[1:6]
  return nothing
end
