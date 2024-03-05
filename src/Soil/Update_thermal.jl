# Campbell 1974, Bonan 2019 Table 8.2
@fastmath function cal_ψ(θ::T, θ_sat::T, ψ_sat::T, b::T) where {T<:Real}
  ψ = ψ_sat * (θ / θ_sat)^(-b)
  max(ψ, ψ_sat)
end

@fastmath cal_K(θ::T, θ_sat::T, K_sat::T, b::T) where {T<:Real} =
  K_sat * (θ / θ_sat)^(2 * b + 3)


# Function to update soil heat flux
function UpdateHeatFlux(p::Soil,
  # Xcs_g::Float64, lambda_snow::Float64, Tsn0::Float64,
  Tair_annual_mean::Float64, period_in_seconds::Float64)

  n = p.n_layer
  # TODO: i may have bug
  @inbounds for i in 2:n+1
    if i <= n
      p.G[i] = (p.Tsoil_p[i-1] - p.Tsoil_p[i]) /
               (0.5 * p.dz[i-1] / p.lambda[i-1] + 0.5 * p.dz[i] / p.lambda[i])
    else
      p.G[i] = p.lambda[i-1] * (p.Tsoil_p[i-1] - Tair_annual_mean) / (DEPTH_F + p.dz[i-1] * 0.5)
    end

    p.G[i] = clamp(p.G[i], -200, 200)
  end

  S = 0.0
  for i in 1:n
    p.Tsoil_c[i] = p.Tsoil_p[i] + (p.G[i] - p.G[i+1] + S) / (p.Cs[i] * p.dz[i]) * period_in_seconds
    p.Tsoil_c[i] = clamp(p.Tsoil_c[i], -50.0, 50.0)
  end

  Update_ice_ratio(p)
  p.Tsoil_p[1:n] .= p.Tsoil_c[1:n]
end


# Function to update volume heat capacity
function Update_Cs(p::Soil)
  for i in 1:p.n_layer
    # Chen B. (2007) Ecological Modelling 209, 277-300  (equation 18)
    term1 = 2.0 * 1.0e+3 * p.density_soil[i] / 2.65
    term2 = 1.0e+6 * p.θ[i] * (4.2 * (1 - p.ice_ratio[i]) + 2.09 * p.ice_ratio[i])
    term3 = 2.5 * 1.0e+6 * p.f_org[i]

    p.Cs[i] = term1 + term2 + term3
  end
end


# Function to update the frozen status of each soil
function Update_ice_ratio(p::Soil)
  Lf0 = 3.34 * 100000  # latent heat of fusion (liquid: solid) at 0C
  # 会不会这里出错了

  @inbounds for i in 1:p.n_layer
    # starting to freeze
    if p.Tsoil_p[i] >= 0.0 && p.Tsoil_c[i] < 0.0 && p.ice_ratio[i] < 1.0
      Gsf = (0.0 - p.Tsoil_c[i]) * p.Cs[i] * p.dz[i]
      p.ice_ratio[i] += Gsf / Lf0 / 1000.0 / (p.θ[i] * p.dz[i])
      p.ice_ratio[i] = min(1.0, p.ice_ratio[i])

      p.Tsoil_c[i] = 0
      # starting to melt
    elseif p.Tsoil_p[i] <= 0.0 && p.Tsoil_c[i] > 0.0 && p.ice_ratio[i] > 0.0
      Gsm = (p.Tsoil_c[i] - 0.0) * p.Cs[i] * p.dz[i]
      p.ice_ratio[i] -= Gsm / Lf0 / 1000.0 / (p.θ[i] * p.dz[i])
      p.ice_ratio[i] = max(0.0, p.ice_ratio[i])

      p.Tsoil_c[i] = 0
    end

    p.ice_ratio[i] *= p.θ_prev[i] / p.θ[i]
    p.ice_ratio[i] = min(1.0, p.ice_ratio[i])
  end
end

# Function to update soil thermal conductivity
function UpdateSoilThermalConductivity(p::Soil)
  ki = 2.1  # the thermal conductivity of ice
  kw = 0.61  # the thermal conductivity of water

  @inbounds for i in 1:p.n_layer
    κ_dry = p.thermal_cond[i]^(1 - p.θ_sat[i])  # dry
    tmp2 = ki^(1.2 * p.θ[i] * p.ice_ratio[i])  # ice.  no source for "1.2"
    tmp3 = kw^(p.θ[i] * (1 - p.ice_ratio[i]))  # water
    tmp4 = p.θ[i] / p.θ_sat[i]  # Sr

    p.lambda[i] = (κ_dry * tmp2 * tmp3 - 0.15) * tmp4 + 0.15  # Note: eq. 8. LHE
    p.lambda[i] = max(p.lambda[i], 0.15)  # juweimin05
  end
end
