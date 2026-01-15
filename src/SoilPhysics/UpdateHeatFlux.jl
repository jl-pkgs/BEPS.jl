# Function to update soil heat flux
function UpdateHeatFlux(st::S, Tair_annual_mean::Float64, period_in_seconds::Float64) where {
  S<:Union{SoilState,Soil}}

  (; G, Tsoil_c, Tsoil_p, κ, dz) = st
  n = st.n_layer

  # TODO: i may have bug
  @inbounds for i in 2:n+1
    if i <= n
      G[i] = 2(Tsoil_p[i-1] - Tsoil_p[i]) / (dz[i-1] / κ[i-1] + dz[i] / κ[i])
    else
      G[i] = κ[i-1] * (Tsoil_p[i-1] - Tair_annual_mean) / (DEPTH_F + dz[i-1] * 0.5)
    end
    G[i] = clamp(G[i], -200, 200)
  end

  source = 0.0
  for i in 1:n
    Tsoil_c[i] = Tsoil_p[i] + (G[i] - G[i+1] + source) / (st.Cv[i] * dz[i]) * period_in_seconds
    Tsoil_c[i] = clamp(Tsoil_c[i], -50.0, 50.0)
  end
  Update_ice_ratio(st)
  Tsoil_p[1:n] .= Tsoil_c[1:n]
end


# Function to update volume heat capacity
# Bonan 2019, Table 5.2
function UpdateThermal_Cv(st::S, ps::P) where {S<:Union{SoilState,Soil},P<:Union{BEPSmodel,Soil}}
  (; θ, ice_ratio) = st
  (; ρ_soil, V_SOM) = get_thermal(ps)

  for i in 1:st.n_layer
    # Chen Baozhang. (2007) Ecological Modelling 209, 277-300  (equation 18)
    term1 = 2.0e+6 * ρ_soil[i] / 2650.0 # soil solid, like Quartz
    term2 = 1.0e+6 * θ[i] * (4.2 * (1 - ice_ratio[i]) + 2.09 * ice_ratio[i]) # water and ice
    term3 = 2.5e+6 * V_SOM[i] # soil organic matter, 2.5 [MJ m-3 K-1]
    st.Cv[i] = term1 + term2 + term3 # [MJ m-3 K-1]
  end
end
UpdateThermal_Cv(p::Soil) = UpdateThermal_Cv(p, p)


# Function to update the frozen status of each soil
# 旧版本：兼容 Soil 结构体
function Update_ice_ratio(st::S) where {S<:Union{SoilState,Soil}}
  (; dz, θ, θ_prev, Tsoil_c, Tsoil_p, ice_ratio, Cv) = st
  Lf0 = 3.34 * 100000  # latent heat of fusion (liquid: solid) at 0C
  # 会不会这里出错了
  @inbounds for i in 1:st.n_layer
    # starting to freeze
    if Tsoil_p[i] >= 0.0 && Tsoil_c[i] < 0.0 && ice_ratio[i] < 1.0
      Gsf = (0.0 - Tsoil_c[i]) * Cv[i] * dz[i]
      ice_ratio[i] += Gsf / Lf0 / 1000.0 / (θ[i] * dz[i])
      ice_ratio[i] = min(1.0, ice_ratio[i])

      Tsoil_c[i] = 0.0
      # starting to melt
    elseif Tsoil_p[i] <= 0.0 && Tsoil_c[i] > 0.0 && ice_ratio[i] > 0.0
      Gsm = (Tsoil_c[i] - 0.0) * Cv[i] * dz[i]
      ice_ratio[i] -= Gsm / Lf0 / 1000.0 / (θ[i] * dz[i])
      ice_ratio[i] = max(0.0, ice_ratio[i])

      Tsoil_c[i] = 0.0
    end

    ice_ratio[i] *= θ_prev[i] / θ[i]
    ice_ratio[i] = min(1.0, ice_ratio[i])
  end
end


# Function to update soil thermal conductivity
@fastmath function UpdateThermal_κ(st::S, ps::P) where {
  S<:Union{SoilState,Soil},P<:Union{BEPSmodel,Soil}}
  (; θ, ice_ratio, κ) = st
  (; θ_sat) = get_hydraulic(ps)
  (; κ_dry) = get_thermal(ps)

  ki = 2.1  # the thermal conductivity of ice
  kw = 0.61  # the thermal conductivity of water

  @inbounds for i in 1:st.n_layer
    k_dry_i = κ_dry[i]^(1 - θ_sat[i])  # dry
    tmp2 = ki^(1.2 * θ[i] * ice_ratio[i])  # ice
    tmp3 = kw^(θ[i] * (1 - ice_ratio[i]))  # water
    tmp4 = θ[i] / θ_sat[i]  # Sr

    κ[i] = (k_dry_i * tmp2 * tmp3 - 0.15) * tmp4 + 0.15  # Note: eq. 8. LHE
    κ[i] = max(κ[i], 0.15) # juweimin05
  end
end
UpdateThermal_κ(p::Soil) = UpdateThermal_κ(p, p)
