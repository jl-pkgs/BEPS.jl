# Function to update soil heat flux
function UpdateHeatFlux(p::Soil,
  # Xcs_g::Float64, κ_snow::Float64, Tsn0::Float64,
  Tair_annual_mean::Float64, period_in_seconds::Float64)
  (; G, Tsoil_c, Tsoil_p, dz, κ) = p

  n = p.n_layer
  # TODO: i may have bug
  @inbounds for i in 2:n+1
    if i <= n
      G[i] = 2(Tsoil_p[i-1] - Tsoil_p[i]) / (dz[i-1] / κ[i-1] + dz[i] / κ[i])
    else
      G[i] = κ[i-1] * (Tsoil_p[i-1] - Tair_annual_mean) / (DEPTH_F + dz[i-1] * 0.5)
    end
    G[i] = clamp(G[i], -200, 200)
  end

  S = 0.0
  for i in 1:n
    Tsoil_c[i] = Tsoil_p[i] + (G[i] - G[i+1] + S) / (p.Cs[i] * dz[i]) * period_in_seconds
    Tsoil_c[i] = clamp(Tsoil_c[i], -50.0, 50.0)
  end

  Update_ice_ratio(p)
  Tsoil_p[1:n] .= Tsoil_c[1:n]
end


# Function to update volume heat capacity
# Bonan 2019, Table 5.2
function Update_Cs(p::Soil)
  (; θ, ice_ratio) = p
  for i in 1:p.n_layer
    # Chen B. (2007) Ecological Modelling 209, 277-300  (equation 18)
    term1 = 2.0 * 1.0e+3 * p.ρ_soil[i] / 2.65 # soil solid
    term2 = 1.0e+6 * θ[i] * (4.2 * (1 - ice_ratio[i]) + 2.09 * ice_ratio[i])
    term3 = 2.5 * 1.0e+6 * p.V_SOM[i] # soil organic matter, 2.5 [MJ m-3 K-1]

    p.Cs[i] = term1 + term2 + term3 # [MJ m-3 K-1]
  end
end


# Function to update the frozen status of each soil
function Update_ice_ratio(p::Soil)
  (; θ, θ_prev, Tsoil_c, Tsoil_p, ice_ratio, dz) = p
  Lf0 = 3.34 * 100000  # latent heat of fusion (liquid: solid) at 0C
  # 会不会这里出错了
  @inbounds for i in 1:p.n_layer
    # starting to freeze
    if Tsoil_p[i] >= 0.0 && Tsoil_c[i] < 0.0 && ice_ratio[i] < 1.0
      Gsf = (0.0 - Tsoil_c[i]) * p.Cs[i] * dz[i]
      ice_ratio[i] += Gsf / Lf0 / 1000.0 / (θ[i] * dz[i])
      ice_ratio[i] = min(1.0, ice_ratio[i])

      Tsoil_c[i] = 0
      # starting to melt
    elseif Tsoil_p[i] <= 0.0 && Tsoil_c[i] > 0.0 && ice_ratio[i] > 0.0
      Gsm = (Tsoil_c[i] - 0.0) * p.Cs[i] * dz[i]
      ice_ratio[i] -= Gsm / Lf0 / 1000.0 / (θ[i] * dz[i])
      ice_ratio[i] = max(0.0, ice_ratio[i])

      Tsoil_c[i] = 0
    end

    ice_ratio[i] *= θ_prev[i] / θ[i]
    ice_ratio[i] = min(1.0, ice_ratio[i])
  end
end

# Function to update soil thermal conductivity
function UpdateSoilThermalConductivity(p::Soil)
  (; θ, ice_ratio, κ, θ_sat) = p
  ki = 2.1  # the thermal conductivity of ice
  kw = 0.61  # the thermal conductivity of water

  @inbounds for i in 1:p.n_layer
    κ_dry = p.κ_dry[i]^(1 - θ_sat[i])  # dry
    tmp2 = ki^(1.2 * θ[i] * ice_ratio[i])  # ice.  no source for "1.2"
    tmp3 = kw^(θ[i] * (1 - ice_ratio[i]))  # water
    tmp4 = θ[i] / θ_sat[i]  # Sr

    κ[i] = (κ_dry * tmp2 * tmp3 - 0.15) * tmp4 + 0.15  # Note: eq. 8. LHE
    κ[i] = max(κ[i], 0.15)  # juweimin05
  end
end

