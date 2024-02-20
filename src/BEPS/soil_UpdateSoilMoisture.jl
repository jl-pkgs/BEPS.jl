# Campbell 1974, Bonan 2019 Table 8.2
function cal_ψ(θ::T, θ_sat::T, ψ_sat::T, b::T) where {T<:Real} 
  ψ = ψ_sat * (θ / θ_sat)^(-b)
  max(ψ, ψ_sat)
end

cal_K(θ::T, θ_sat::T, K_sat::T, b::T) where {T<:Real} = K_sat * (θ / θ_sat)^(2 * b + 3)


function UpdateSoilMoisture(p::Soil, kstep::Float64)
  inf, inf_max = 0.0, 0.0
  this_step, total_t, max_Fb = 0.0, 0.0, 0.0

  @unpack f_ice, Ksat, Kb, KK, km, b = p
  dz = p.d_soil
  θ = p.thetam
  θb = p.thetab
  n = p.n_layer
  ψ_sat = p.psi_sat
  ψ = p.psim
  ψb = p.psib
  θ_sat = p.fei # porosity
  p.thetam_prev .= p.thetam # assign current soil moisture to prev
  θ_prev = p.thetam_prev

  # TODO: check this
  @inbounds for i in 1:n+1
    if p.temp_soil_c[i] > 0.0
      f_ice[i] = 1.0
    elseif p.temp_soil_c[i] < -1.0
      f_ice[i] = 0.1
    else
      f_ice[i] = 0.1 + 0.9 * (p.temp_soil_c[i] + 1.0)
    end
  end

  # Max infiltration calculation
  inf_max = f_ice[1] * Ksat[1] * (1 + (θ_sat[1] - θ_prev[1]) / dz[1] * ψ_sat[1] * b[1] / θ_sat[1])
  inf = max(f_ice[1] * (p.Zp / kstep + p.r_rain_g), 0)
  inf = clamp(inf, 0, inf_max)

  # Ponded water after runoff. This one is related to runoff. LHe.
  p.Zp = (p.Zp / kstep + p.r_rain_g - inf) * kstep * p.r_drainage

  @inbounds while total_t < kstep

    # soil moisture in the boundaries
    for i in 1:n
      if i < n
        θb[i] = (θ[i+1] / dz[i+1] + θ[i] / dz[i]) / (1 / dz[i] + 1 / dz[i+1])
      else
        d1 = max((θ[i] - θb[i-1]) * 2.0 / dz[i], 0)
        θb[i] = min(θ[i] + d1 * dz[i] / 2.0, θ_sat[i])
      end
    end

    for i in 1:n-1
      K = f_ice[i] * (θb[i] / θ_sat[i])^(2 * b[i] + 3)
      Kb[i] = (Ksat[i] * dz[i] + Ksat[i+1] * dz[i+1]) / (dz[i] + dz[i+1]) * K
    end
    Kb[n] = 0.5 * f_ice[n] * cal_K(θ[n], θ_sat[n], Ksat[n], b[n])

    # the unsaturated soil water retention. LHe
    for i in 1:n
      ψ[i] = cal_ψ(θ[i], θ_sat[i], ψ_sat[i], b[i])
      ψb[i] = cal_ψ(θb[i], θ_sat[i], ψ_sat[i], b[i])
    end
    
    # Hydraulic conductivity: Bonan, Table 8.2, Campbell 1974, K = K_sat*(θ/θ_sat)^(2b+3)
    for i in 1:n
      km[i] = f_ice[i] * cal_K(θ[i], θ_sat[i], Ksat[i], b[i])
    end

    # the unsaturated hydraulic conductivity of soil layer @ boundaries
    for i in 1:n-1
      KK[i] = (km[i] * ψ[i] + km[i+1] * ψ[i+1]) / (ψ[i] + ψ[i+1]) * (b[i] + b[i+1]) / (b[i] + b[i+1] + 6)
    end
    KK[n] = (km[n] * ψ[n] + Kb[n] * ψb[n]) / (ψ[n] + ψb[n]) * b[n] / (b[n] + 3)

    # Fb, flow speed. Dancy's law. LHE.
    for i in 1:n-1
      p.r_waterflow[i] = KK[i] * (2 * (ψ[i+1] - ψ[i]) / (dz[i] + dz[i+1]) + 1)
    end
    p.r_waterflow[n] = 0

    # check the r_waterflow further. LHE
    for i in 1:n-1
      p.r_waterflow[i] = min((θ_sat[i+1] - θ[i+1]) * dz[i+1] / kstep + p.Ett[i+1], p.r_waterflow[i])
      max_Fb = max(max_Fb, abs(p.r_waterflow[i]))
    end

    this_step = guess_step(max_Fb)
    total_t += this_step
    total_t > kstep && (this_step -= (total_t - kstep))

    # from there: kstep is replaced by this_step. LHE
    for i in 1:n
      if i == 1
        θ[i] += (inf - p.r_waterflow[i] - p.Ett[i]) * this_step / dz[i]
      else
        θ[i] += (p.r_waterflow[i-1] - p.r_waterflow[i] - p.Ett[i]) * this_step / dz[i]
      end
      θ[i] = clamp(θ[i], p.theta_vwp[i], θ_sat[i])
    end
  end

  for i in 1:n
    p.ice_ratio[i] *= θ_prev[i] / θ[i]
    p.ice_ratio[i] = min(1.0, p.ice_ratio[i])
  end
end



function guess_step(max_Fb)
  if max_Fb > 1.0e-5
    this_step = 1.0
  elseif max_Fb > 1.0e-6
    this_step = 30.0
  else
    this_step = 360.0
  end
  this_step
end
