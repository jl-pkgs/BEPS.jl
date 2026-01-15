export Init_Soil_Parameters, Init_Soil_T_θ!, Sync_Soil_to_State!, UpdateRootFraction, Update_Tsoil_c, Update_G
export UpdateHeatFlux, UpdateThermal_Cs,
  Update_ice_ratio,
  UpdateThermal_κ,
  soil_water_factor_v2,
  UpdateSoilMoisture, Root_Water_Uptake


get_hydraulic(ps::BEPSmodel) = ps.hydraulic
get_hydraulic(ps::Soil) = ps

get_thermal(ps::BEPSmodel) = ps.thermal
get_thermal(ps::Soil) = ps

get_root_decay(ps::BEPSmodel) = ps.veg.r_root_decay
get_root_decay(ps::Soil) = ps.r_root_decay

include("UpdateHeatFlux.jl")
include("UpdateSoilMoisture.jl")
include("soil_water_factor_v2.jl")


function Update_Tsoil_c(st::S, value::Cdouble) where {S<:Union{SoilState,Soil}}
  st.Tsoil_c[1] = value
end

function Update_G(st::S, value::Cdouble) where {S<:Union{SoilState,Soil}}
  st.G[1] = value
end

# T < -1℃, all frozen; T > 0℃, no frozen; else partially frozen
get_ice_ratio(Tsoil::FT) where {FT} = clamp(-Tsoil, FT(0), FT(1))

# 新版本：JAX 风格 (st, ps) 签名
function UpdateRootFraction!(st::S, ps::P) where {
  S<:Union{SoilState,Soil},P<:Union{BEPSmodel,Soil}}

  n = st.n_layer
  (; f_root, dz) = st
  β = get_root_decay(ps)

  z = zeros(n) # cumulative depth of soil layers
  z[1] = dz[1] * 100
  f_root[1] = 1 - β^(z[1])

  for i in 2:(n-1)
    z[i] = z[i-1] + dz[i] * 100
    f_root[i] = β^(z[i-1]) - β^(z[i])
  end
  f_root[n] = β^(z[n-1])
end
UpdateRootFraction!(soil::Soil) = UpdateRootFraction!(soil, soil)


"""
    Init_Soil_T_θ!(p::Soil, Tsoil, Tair, θ0, snowdepth)

初始化土壤剖面的温度(T)和含水量(θ)。
"""
function Init_Soil_T_θ!(st::S, Tsoil::Float64, Tair::Float64, θ0::Float64, snowdepth::Float64) where {
  S<:Union{SoilState,Soil}}

  dT = clamp(Tsoil - Tair, -5.0, 5.0)
  st.z_snow = snowdepth
  # p.z_water = 0.0
  # p.r_rain_g = 0.0
  T_scale_factors = [0.4, 0.5, 1.0, 1.2, 1.4]
  θ_scale_factors = [0.8, 1.0, 1.05, 1.10, 1.15]

  for i in 1:st.n_layer
    st.Tsoil_c[i] = Tair + T_scale_factors[i] * dT
    st.Tsoil_p[i] = Tair + T_scale_factors[i] * dT
    st.θ[i] = θ_scale_factors[i] * θ0
    st.θ_prev[i] = θ_scale_factors[i] * θ0
    st.ice_ratio[i] = get_ice_ratio(st.Tsoil_c[i])
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
