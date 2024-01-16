# Function to update soil heat flux
function UpdateHeatFlux(p::Soil,
  Xg_snow::Float64, lambda_snow::Float64, Tsn0::Float64,
  Tair_annual_mean::Float64, period_in_seconds::Float64)

  # TODO: i may have bug
  @inbounds for i in 2:p.n_layer
    if i <= p.n_layer
      p.G[i] = (p.temp_soil_p[i-1] - p.temp_soil_p[i]) / (0.5 * p.d_soil[i-1] / p.lambda[i-1] + 0.5 * p.d_soil[i] / p.lambda[i])
    else
      p.G[i] = p.lambda[i-1] * (p.temp_soil_p[i-1] - Tair_annual_mean) / (DEPTH_F + p.d_soil[i-1] * 0.5)
    end

    p.G[i] = clamp(p.G[i], -200, 200)
  end

  S = 0.0
  for i in 1:p.n_layer
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
  ft = zeros(Float64, p.n_layer)
  fpsisr = zeros(Float64, p.n_layer)
  dtt = zeros(Float64, p.n_layer)

  t1 = -0.02
  t2 = 2.0

  if p.psim[1] <= 0.000001
    for i in 1:p.n_layer
      p.psim[i] = p.psi_sat[i] * (p.thetam[i] / p.fei[i])^-p.b[i]
      p.psim[i] = max(p.psi_sat[i], p.psim[i])
    end
  end

  for i in 1:p.n_layer
    # psi_sr in m H2O! This is the old version. LHE.
    fpsisr[i] = p.psim[i] > p.psi_min ? 1.0 / (1 + ((p.psim[i] - p.psi_min) / p.psi_min)^p.alpha) : 0

    ft[i] = p.temp_soil_p[i] > 0.0 ? 1.0 - exp(t1 * p.temp_soil_p[i]^t2) : 0

    fpsisr[i] *= ft[i]
  end

  for i in 1:p.n_layer
    dtt[i] = FW_VERSION == 1 ? p.f_root[i] * fpsisr[i] : p.f_root[i]
  end
  dtt_sum = sum(dtt)

  if dtt_sum < 0.000001
    p.f_soilwater = 0.1
  else
    for i in 1:p.n_layer
      p.dt[i] = dtt[i] / dtt_sum

      if isnan(p.dt[i])
        println(p.dt[1])
      end
    end

    fpsisr_sum = sum(fpsisr .* p.dt)
    p.f_soilwater = max(0.1, fpsisr_sum)
  end
end


function UpdateSoilMoisture(p::Soil, kstep::Float64)
  inf, inf_max = 0.0, 0.0
  this_step, total_t, max_Fb = 0.0, 0.0, 0.0

  # assign the current soil temperature to prev variables.
  p.thetam_prev .= p.thetam

  # Compute f_ice
  @inbounds for i in 1:p.n_layer
    if p.temp_soil_c[i] > 0.0
      p.f_ice[i] = 1.0
    elseif p.temp_soil_c[i] < -1.0
      p.f_ice[i] = 0.1
    else
      p.f_ice[i] = 0.1 + 0.9 * (p.temp_soil_c[i] + 1.0)
    end
  end

  # Max infiltration calculation
  inf_max = p.f_ice[1] * p.Ksat[1] * (1 + (p.fei[1] - p.thetam_prev[1]) / p.d_soil[1] * p.psi_sat[1] * p.b[1] / p.fei[1])
  inf = max(p.f_ice[1] * (p.Zp / kstep + p.r_rain_g), 0)
  inf = clamp(inf, 0, inf_max)

  # Ponded water after runoff. This one is related to runoff. LHe.
  p.Zp = (p.Zp / kstep + p.r_rain_g - inf) * kstep * p.r_drainage

  @inbounds while total_t < kstep
    for i in 1:p.n_layer
      p.km[i] = p.f_ice[i] * p.Ksat[i] * ((p.thetam[i] / p.fei[i])^(2 * p.b[i] + 3))
    end

    # soil moisture in the boundaries
    for i in 1:p.n_layer
      if i < p.n_layer
        p.thetab[i] = (p.thetam[i+1] / p.d_soil[i+1] + p.thetam[i] / p.d_soil[i]) / (1 / p.d_soil[i] + 1 / p.d_soil[i+1])
      else
        d1 = max((p.thetam[i] - p.thetab[i-1]) * 2.0 / p.d_soil[i], 0)
        p.thetab[i] = p.thetam[i] + d1 * p.d_soil[i] / 2.0
        p.thetab[i] = min(p.thetab[i], p.fei[i])
      end
    end

    for i in 1:p.n_layer
      if i < p.n_layer
        p.Kb[i] = p.f_ice[i] * (p.Ksat[i] * p.d_soil[i] + p.Ksat[i+1] * p.d_soil[i+1]) / (p.d_soil[i] + p.d_soil[i+1]) * (p.thetab[i] / p.fei[i])^(2 * p.b[i] + 3)
      else
        p.Kb[i] = 0.5 * p.f_ice[i] * p.Ksat[i] * (p.thetab[i] / p.fei[i])^(2 * p.b[i] + 3)
      end
    end

    # the unsaturated soil water retention. LHe
    for i in 1:p.n_layer
      p.psim[i] = p.psi_sat[i] * (p.thetam[i] / p.fei[i])^(-p.b[i])
      p.psim[i] = max(p.psi_sat[i], p.psim[i])
    end

    # the unsaturated soil water retention @ boundaries. LHe
    for i in 1:p.n_layer
      p.psib[i] = p.psi_sat[i] * (p.thetab[i] / p.fei[i])^(-p.b[i])
      p.psib[i] = max(p.psi_sat[i], p.psib[i])
    end

    # the unsaturated hydraulic conductivity of soil layer @ boundaries
    for i in 1:p.n_layer
      if i < p.n_layer
        p.KK[i] = (p.km[i] * p.psim[i] + p.km[i+1] * p.psim[i+1]) / (p.psim[i] + p.psim[i+1]) * (p.b[i] + p.b[i+1]) / (p.b[i] + p.b[i+1] + 6)
      else
        p.KK[i] = (p.km[i] * p.psim[i] + p.Kb[i] * p.psib[i]) / (p.psim[i] + p.psib[i]) * p.b[i] / (p.b[i] + 3)
      end
    end

    # Fb, flow speed. Dancy's law. LHE.
    for i in 1:p.n_layer
      if i < p.n_layer
        p.r_waterflow[i] = p.KK[i] * (2 * (p.psim[i+1] - p.psim[i]) / (p.d_soil[i] + p.d_soil[i+1]) + 1)
      else
        p.r_waterflow[i] = 0
      end
    end

    # check the r_waterflow further. LHE
    for i in 1:p.n_layer-1
      p.r_waterflow[i] = min((p.fei[i+1] - p.thetam[i+1]) * p.d_soil[i+1] / kstep + p.Ett[i+1], p.r_waterflow[i])
      max_Fb = max(max_Fb, abs(p.r_waterflow[i]))
    end

    if max_Fb > 1.0e-5
      this_step = 1.0
    elseif max_Fb > 1.0e-6
      this_step = 30.0
    else
      this_step = 360.0
    end

    total_t += this_step
    total_t > kstep && (this_step -= (total_t - kstep))

    # from there: kstep is replaced by this_step. LHE
    for i in 1:p.n_layer
      if i == 1
        p.thetam[i] += (inf * this_step - p.r_waterflow[i] * this_step - p.Ett[i] * this_step) / p.d_soil[i]
      else
        p.thetam[i] += (p.r_waterflow[i-1] * this_step - p.r_waterflow[i] * this_step - p.Ett[i] * this_step) / p.d_soil[i]
      end

      p.thetam[i] = clamp(p.thetam[i], p.theta_vwp[i], p.fei[i])
    end
  end

  for i in 1:p.n_layer
    p.ice_ratio[i] *= p.thetam_prev[i] / p.thetam[i]
    p.ice_ratio[i] = min(1.0, p.ice_ratio[i])
  end
end


# Function to calculate soil water uptake from a layer
function Soil_Water_Uptake(p::Soil, Trans_o::Float64, Trans_u::Float64, Evap_soil::Float64)
  rho_w = 1025.0
  Trans = Trans_o + Trans_u

  # for the top layer
  p.Ett[1] = Trans / rho_w * p.dt[1] + Evap_soil / rho_w

  # for each layer:
  for i in 2:p.n_layer
    p.Ett[i] = Trans / rho_w * p.dt[i]
  end
end

export UpdateHeatFlux, Update_Cs, UpdateSoilThermalConductivity, soil_water_factor_v2, UpdateSoilMoisture, Soil_Water_Uptake
