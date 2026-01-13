# landcover::Int, soil_type::Int, Tsoil, soilwater, snowdepth
# initialize soil conditions, read soil parameters and set depth
function Init_Soil_var(soil::AbstractSoil, state::Union{State,Vector}, Ta::FT;
  VegType::Int=25, SoilType::Int=8,
  r_drainage::FT, r_root_decay::FT,
  Tsoil0::FT=Ta, θ0::FT=0.2, z_snow0::FT=0.0,
  ignored...
) where {FT<:Real}
  # r_drainage = param[27]
  # r_root_decay = param[28]
  soil.r_drainage = r_drainage
  Init_Soil_Parameters(VegType, SoilType, r_root_decay, soil)
  Init_Soil_Status(soil, Tsoil0, Ta, θ0, z_snow0) # LHE
  Init_State!(state, soil, Ta)
end

# for (i=3;i<=8;i++)   var_o[i] = tem;
# for(i=9;i<=14;i++) var_o[i] = soil->Tsoil_p[i-9];
# for(i=21;i<=26;i++) var_o[i] = soil->θ_prev[i-21];
# for(i=27;i<=32;i++) var_o[i] = soil->ice_ratio[i-27];
function Init_State!(state::Vector, soil::AbstractSoil, Ta)
  state .= 0
  for i = 1:6
    state[i+3] = Ta
    state[i+9] = soil.Tsoil_p[i]
    state[i+21] = soil.θ_prev[i]
    state[i+27] = soil.ice_ratio[i]
  end
  return nothing
end

function Init_State!(state::State, soil::AbstractSoil, Ta)
  state.Ts .= Ta
  state.Tsoil_prev .= soil.Tsoil_p[1:6]
  state.θ_prev .= soil.θ_prev[1:6]
  state.ice_ratio .= soil.ice_ratio[1:6]
  return nothing
end

# function update_state!(state::State, var_n::Vector{Float64})
#   Tsoil = state.Tsoil
#   var_n[3+1] = Tsoil[1]       # Ts0, 3
#   var_n[4+1] = Tsoil[2]       # Tsn, 4 
#   var_n[5+1] = Tsoil[3]       # Tsm0, 5
#   var_n[6+1] = Tsoil[4]       # Tsn1, 6
#   var_n[7+1] = Tsoil[5]       # Tsn1, 7

#   var_n[11+1] = state.Qhc_o   # Qhc_o, 11, sensible heat flux

#   var_n[15+1] = state.m_water.o
#   var_n[18+1] = state.m_water.u

#   var_n[16+1] = state.m_snow.o
#   var_n[19+1] = state.m_snow.u
#   var_n[20+1] = state.m_snow.g
# end
