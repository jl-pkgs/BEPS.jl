# T < -1℃, all frozen; T > 0℃, no frozen; else partially frozen
get_ice_ratio(Tsoil::FT) where {FT} = clamp(-Tsoil, FT(0), FT(1))


# update `f_root`
function UpdateRootFraction!(soil::Soil)
  f_root = soil.f_root
  β = soil.r_root_decay

  cum_depth = zeros(soil.n_layer) # cumulative depth of soil layers
  cum_depth[1] = soil.dz[1] 
  f_root[1] = 1 - β^(cum_depth[1] * 100)

  # For 2 to n_layer - 1
  for i in 2:(soil.n_layer-1)
    cum_depth[i] = cum_depth[i-1] + soil.dz[i]
    f_root[i] = β^(cum_depth[i-1] * 100) - β^(cum_depth[i] * 100)
  end
  # For the last layer
  f_root[soil.n_layer] = β^(cum_depth[soil.n_layer-1] * 100)
end


"""
    Init_Soil_T_θ!(p::Soil, Tsoil, Tair, θ0, snowdepth)

初始化土壤剖面的温度(T)和含水量(θ)。
"""
function Init_Soil_T_θ!(p::Soil, Tsoil::Float64, Tair::Float64, θ0::Float64, snowdepth::Float64)
  d_t = clamp(Tsoil - Tair, -5.0, 5.0)
  # p.z_water = 0.0
  # p.r_rain_g = 0.0
  p.z_snow = snowdepth

  temp_scale_factors = [0.4, 0.5, 1.0, 1.2, 1.4]
  moisture_scale_factors = [0.8, 1.0, 1.05, 1.10, 1.15]

  for i in 1:p.n_layer
    p.Tsoil_c[i] = Tair + temp_scale_factors[i] * d_t
    p.Tsoil_p[i] = Tair + temp_scale_factors[i] * d_t
    p.θ[i] = moisture_scale_factors[i] * θ0
    p.θ_prev[i] = moisture_scale_factors[i] * θ0
    p.ice_ratio[i] = get_ice_ratio(p.Tsoil_c[i])
  end
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
