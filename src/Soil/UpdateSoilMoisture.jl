# LSM of Xuanze Zhang
# 
# Soil moisture is predicted from a 5-layer model (as with soil
# temperature), in which the vertical soil moisture transport is governed
# by infiltration, runoff, gradient diffusion, gravity, and root
# extraction through canopy transpiration.  The net water applied to the
# surface layer is the snowmelt plus precipitation plus the throughfall
# of canopy dew minus surface runoff and evaporation.
# CLM3.5 uses a zero-flow bottom boundary condition.
function UpdateSoilMoisture(soil::Soil, kstep::Float64)
  inf, inf_max = 0.0, 0.0
  Δt, total_t, max_Fb = 0.0, 0.0, 0.0

  n = soil.n_layer
  @unpack dz, f_water, Ksat, KK, km, b,
  ψ_sat, ψ,
  θ_sat, θ, θ_prev = soil
  θ_prev .= θ # assign current soil moisture to prev

  # ψb, θb,
  # TODO: check this
  @inbounds for i in 1:n+1
    if soil.Tsoil_c[i] > 0.0
      f_water[i] = 1.0
    elseif soil.Tsoil_c[i] < -1.0
      f_water[i] = 0.1
    else
      f_water[i] = 0.1 + 0.9 * (soil.Tsoil_c[i] + 1.0) # 冰含量的一个指标
    end
  end

  # Max infiltration calculation
  # Ksat * (1 + (θ_sat - θ_prev) / dz * ψ_sat / θ_sat * b )
  inf_max = f_water[1] * Ksat[1] * (1 + (θ_sat[1] - θ_prev[1]) / dz[1] * ψ_sat[1] * b[1] / θ_sat[1])
  inf = max(f_water[1] * (soil.z_water / kstep + soil.r_rain_g), 0)
  inf = clamp(inf, 0, inf_max)

  # Ponded water after runoff. This one is related to runoff. LHe.
  soil.z_water = (soil.z_water / kstep + soil.r_rain_g - inf) * kstep * soil.r_drainage

  @inbounds while total_t < kstep
    # 为了解决相互依赖的关系，循环寻找稳态
    # the unsaturated soil water retention. LHe
    # Hydraulic conductivity: Bonan, Table 8.2, Campbell 1974, K = K_sat*(θ/θ_sat)^(2b+3)
    for i in 1:n
      ψ[i] = cal_ψ(θ[i], θ_sat[i], ψ_sat[i], b[i])
      km[i] = f_water[i] * cal_K(θ[i], θ_sat[i], Ksat[i], b[i]) # Hydraulic conductivity, [m/s]
    end

    # Fb, flow speed. Dancy's law. LHE.
    # check the r_waterflow further. LHE
    for i in 1:n-1
      # 不同层土壤深度不同，能否这样写？
      # K * ψ * b / (b + 3): ?
      # the unsaturated hydraulic conductivity of soil layer
      KK[i] = (km[i] * ψ[i] + km[i+1] * ψ[i+1]) / (ψ[i] + ψ[i+1]) * (b[i] + b[i+1]) / (b[i] + b[i+1] + 6) # 计算平均的一种方案？
      Q = KK[i] * (2 * (ψ[i+1] - ψ[i]) / (dz[i] + dz[i+1]) + 1) # z direction
      Q_max = (θ_sat[i+1] - θ[i+1]) * dz[i+1] / kstep + soil.Ett[i+1]
      Q = min(Q, Q_max)

      soil.r_waterflow[i] = Q
      max_Fb = max(max_Fb, abs(Q))
    end
    # p.r_waterflow[n] = 0

    Δt = guess_step(max_Fb) # this_step
    total_t += Δt
    total_t > kstep && (Δt -= (total_t - kstep))

    # from there: kstep is replaced by this_step. LHE
    for i in 1:n
      if i == 1
        θ[i] += (inf - soil.r_waterflow[i] - soil.Ett[i]) * Δt / dz[i]
      else
        θ[i] += (soil.r_waterflow[i-1] - soil.r_waterflow[i] - soil.Ett[i]) * Δt / dz[i]
      end
      θ[i] = clamp(θ[i], soil.θ_vwp[i], θ_sat[i])
    end
  end

  for i in 1:n
    soil.ice_ratio[i] *= θ_prev[i] / θ[i]
    soil.ice_ratio[i] = min(1.0, soil.ice_ratio[i])
  end
end


# Campbell 1974, Bonan 2019 Table 8.2
@fastmath function cal_ψ(θ::T, θ_sat::T, ψ_sat::T, b::T) where {T<:Real}
  ψ = ψ_sat * (θ / θ_sat)^(-b)
  max(ψ, ψ_sat)
end

@fastmath cal_K(θ::T, θ_sat::T, K_sat::T, b::T) where {T<:Real} =
  K_sat * (θ / θ_sat)^(2 * b + 3)


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

# Function to calculate soil water uptake from a layer
"""
    Root Water Uptake

- `土壤蒸发`：仅发生在表层
- `植被蒸腾`：根据根系分布，耗水可能来自于土壤的每一层
"""
function Root_Water_Uptake(p::Soil, Trans_o::Float64, Trans_u::Float64, Evap_soil::Float64)
  Trans = Trans_o + Trans_u
  p.Ett[1] = Trans / ρ_w * p.dt[1] + Evap_soil / ρ_w # for the top layer
  for i in 2:p.n_layer
    p.Ett[i] = Trans / ρ_w * p.dt[i]
  end
end
