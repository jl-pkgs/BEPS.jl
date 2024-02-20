# Function to update soil heat flux
function UpdateHeatFlux(p::Soil,
  Xg_snow::Float64, lambda_snow::Float64, Tsn0::Float64,
  Tair_annual_mean::Float64, period_in_seconds::Float64)

  n = p.n_layer
  # TODO: i may have bug
  @inbounds for i in 2:n+1
    if i <= n
      p.G[i] = (p.temp_soil_p[i-1] - p.temp_soil_p[i]) / (0.5 * p.d_soil[i-1] / p.lambda[i-1] + 0.5 * p.d_soil[i] / p.lambda[i])
    else
      p.G[i] = p.lambda[i-1] * (p.temp_soil_p[i-1] - Tair_annual_mean) / (DEPTH_F + p.d_soil[i-1] * 0.5)
    end

    p.G[i] = clamp(p.G[i], -200, 200)
  end

  S = 0.0
  for i in 1:n
    p.temp_soil_c[i] = p.temp_soil_p[i] + (p.G[i] - p.G[i+1] + S) / (p.Cs[i] * p.d_soil[i]) * period_in_seconds
    p.temp_soil_c[i] = clamp(p.temp_soil_c[i], -50.0, 50.0)
  end

  Update_ice_ratio(p)
  p.temp_soil_p .= p.temp_soil_c
end


# Function to update volume heat capacity
function Update_Cs(p::Soil)
  for i in 1:p.n_layer
    # Chen B. (2007) Ecological Modelling 209, 277-300  (equation 18)
    term1 = 2.0 * 1.0e+3 * p.density_soil[i] / 2.65
    term2 = 1.0e+6 * p.thetam[i] * (4.2 * (1 - p.ice_ratio[i]) + 2.09 * p.ice_ratio[i])
    term3 = 2.5 * 1.0e+6 * p.f_org[i]

    p.Cs[i] = term1 + term2 + term3
  end
end


# Function to update the frozen status of each soil
function Update_ice_ratio(p::Soil)
  Lf0 = 3.34 * 100000  # latent heat of fusion (liquid: solid) at 0C

  @inbounds for i in 1:p.n_layer
    # starting to freeze
    if p.temp_soil_p[i] >= 0.0 && p.temp_soil_c[i] < 0.0 && p.ice_ratio[i] < 1.0
      Gsf = (0.0 - p.temp_soil_c[i]) * p.Cs[i] * p.d_soil[i]
      p.ice_ratio[i] += Gsf / Lf0 / 1000.0 / (p.thetam[i] * p.d_soil[i])
      p.ice_ratio[i] = min(1.0, p.ice_ratio[i])

      p.temp_soil_c[i] = 0
      # starting to melt
    elseif p.temp_soil_p[i] <= 0.0 && p.temp_soil_c[i] > 0.0 && p.ice_ratio[i] > 0.0
      Gsm = (p.temp_soil_c[i] - 0.0) * p.Cs[i] * p.d_soil[i]
      p.ice_ratio[i] -= Gsm / Lf0 / 1000.0 / (p.thetam[i] * p.d_soil[i])
      p.ice_ratio[i] = max(0.0, p.ice_ratio[i])

      p.temp_soil_c[i] = 0
    end

    p.ice_ratio[i] *= p.thetam_prev[i] / p.thetam[i]
    p.ice_ratio[i] = min(1.0, p.ice_ratio[i])
  end
end

# Function to update soil thermal conductivity
function UpdateSoilThermalConductivity(p::Soil)
  ki = 2.1  # the thermal conductivity of ice
  kw = 0.61  # the thermal conductivity of water

  @inbounds for i in 1:p.n_layer
    tmp1 = p.thermal_cond[i]^(1 - p.fei[i])  # dry
    tmp2 = ki^(1.2 * p.thetam[i] * p.ice_ratio[i])  # ice.  no source for "1.2"
    tmp3 = kw^(p.thetam[i] * (1 - p.ice_ratio[i]))  # water
    tmp4 = p.thetam[i] / p.fei[i]  # Sr

    p.lambda[i] = (tmp1 * tmp2 * tmp3 - 0.15) * tmp4 + 0.15  # Note: eq. 8. LHE
    p.lambda[i] = max(p.lambda[i], 0.15)  # juweimin05
  end
end


# Function to compute soil water stress factor
function soil_water_factor_v2(p::Soil)
  # @unpack psim, psi_sat, 
  θ = p.thetam
  n = p.n_layer

  t1 = -0.02
  t2 = 2.0

  if p.psim[1] <= 0.000001
    for i in 1:n
      p.psim[i] = p.psi_sat[i] * (θ[i] / p.fei[i])^-p.b[i]
      p.psim[i] = max(p.psi_sat[i], p.psim[i])
    end
  end

  for i in 1:n
    # psi_sr in m H2O! This is the old version. LHE.
    p.fpsisr[i] = p.psim[i] > p.psi_min ? 1.0 / (1 + ((p.psim[i] - p.psi_min) / p.psi_min)^p.alpha) : 1.0

    p.ft[i] = p.temp_soil_p[i] > 0.0 ? 1.0 - exp(t1 * p.temp_soil_p[i]^t2) : 0

    p.fpsisr[i] *= p.ft[i]
  end

  for i in 1:n
    p.dtt[i] = FW_VERSION == 1 ? p.f_root[i] * p.fpsisr[i] : p.f_root[i]
  end
  dtt_sum = sum(p.dtt)

  if dtt_sum < 0.000001
    p.f_soilwater = 0.1
  else
    fpsisr_sum = 0
    for i in 1:n
      p.dt[i] = p.dtt[i] / dtt_sum

      if isnan(p.dt[i])
        println(p.dt[1])
      end
      fpsisr_sum += p.fpsisr[i] * p.dt[i]
    end
    p.f_soilwater = max(0.1, fpsisr_sum)
  end
end


# Function to calculate soil water uptake from a layer
function Soil_Water_Uptake(p::Soil, Trans_o::Float64, Trans_u::Float64, Evap_soil::Float64)
  ρ_w = 1025.0
  Trans = Trans_o + Trans_u

  # for the top layer
  p.Ett[1] = Trans / ρ_w * p.dt[1] + Evap_soil / ρ_w

  # for each layer:
  for i in 2:p.n_layer
    p.Ett[i] = Trans / ρ_w * p.dt[i]
  end
end

include("soil_UpdateSoilMoisture.jl")

export UpdateHeatFlux, Update_Cs,
  Update_ice_ratio,
  UpdateSoilThermalConductivity,
  soil_water_factor_v2,
  UpdateSoilMoisture,
  Soil_Water_Uptake
