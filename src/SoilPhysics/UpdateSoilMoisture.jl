# LSM of Xuanze Zhang
# 
# Soil moisture is predicted from a 5-layer model (as with soil
# temperature), in which the vertical soil moisture transport is governed
# by infiltration, runoff, gradient diffusion, gravity, and root
# extraction through canopy transpiration.  The net water applied to the
# surface layer is the snowmelt plus precipitation plus the throughfall
# of canopy dew minus surface runoff and evaporation.
# CLM3.5 uses a zero-flow bottom boundary condition.

# 旧版本：兼容 Soil 结构体
UpdateSoilMoisture(soil::Soil, kstep::Float64) = UpdateSoilMoisture(soil, soil, kstep)

# 新版本：JAX 风格 (st, ps) 签名
function UpdateSoilMoisture(st::S, ps::P, kstep::Float64) where {
  S<:Union{StateBEPS,Soil},P<:Union{ParamBEPS,Soil}}

  n = st.n_layer
  dz = st.dz

  r_drainage = ps.r_drainage
  (; θ_sat, K_sat, ψ_sat, b, θ_vwp) = get_hydraulic(ps)
  (; f_water, KK, km, ψ, θ, θ_prev, Tsoil_c, Ett, r_waterflow, ice_ratio, z_water, r_rain_g) = st

  θ_prev .= θ

  @inbounds for i in 1:n+1
    # 注意：Tsoil_c 长度通常是 n，但这里循环到 n+1，需确认 Tsoil_c 实际分配长度。
    # 假设 Tsoil_c 长度足够，或者边界处理
    # Soil struct 定义 dz 为 Vector{Float64} = zeros(10)，所以可以到 n+1 (5+1=6)
    val_T = i <= length(Tsoil_c) ? Tsoil_c[i] : Tsoil_c[end] # 简单边界保护

    if val_T > 0.0
      f_water[i] = 1.0
    elseif val_T < -1.0
      f_water[i] = 0.1
    else
      f_water[i] = 0.1 + 0.9 * (val_T + 1.0)
    end
  end

  # Max infiltration calculation
  # K_sat * (1 + (θ_sat - θ_prev) / dz * ψ_sat / θ_sat * b)
  inf_max = f_water[1] * K_sat[1] * (1 + (θ_sat[1] - θ_prev[1]) / dz[1] * ψ_sat[1] * b[1] / θ_sat[1])
  inf = max(f_water[1] * (z_water / kstep + r_rain_g), 0)
  inf = clamp(inf, 0, inf_max)

  # Ponded water after runoff
  st.z_water = (z_water / kstep + r_rain_g - inf) * kstep * r_drainage

  total_t, max_Fb = 0.0, 0.0

  @inbounds while total_t < kstep
    # 为了解决相互依赖的关系，循环寻找稳态
    # the unsaturated soil water retention. LHe
    # Hydraulic conductivity: Bonan, Table 8.2, Campbell 1974, K = K_sat*(θ/θ_sat)^(2b+3)
    for i in 1:n
      ψ[i] = cal_ψ(θ[i], θ_sat[i], ψ_sat[i], b[i])
      km[i] = f_water[i] * cal_K(θ[i], θ_sat[i], K_sat[i], b[i]) # Hydraulic conductivity, [m/s]
    end

    # Fb, flow speed. Dancy's law. LHE.
    # check the r_waterflow further. LHE
    for i in 1:n-1
      # 不同层土壤深度不同，能否这样写？
      # K * ψ * b / (b + 3): ?
      # the unsaturated hydraulic conductivity of soil layer
      KK[i] = (km[i] * ψ[i] + km[i+1] * ψ[i+1]) / (ψ[i] + ψ[i+1]) * (b[i] + b[i+1]) / (b[i] + b[i+1] + 6) # 计算平均的一种方案？
      Q = KK[i] * (2 * (ψ[i+1] - ψ[i]) / (dz[i] + dz[i+1]) + 1) # z direction
      # `Q_max`出现了单位不匹配的问题，导致Q_max未发挥作用
      Q_max = (θ_sat[i+1] - θ[i+1]) * dz[i+1] / kstep + Ett[i+1]
      Q = min(Q, Q_max)

      r_waterflow[i] = Q
      max_Fb = max(max_Fb, abs(Q))
    end
    # p.r_waterflow[n] = 0

    Δt = guess_step(max_Fb) # this_step
    total_t += Δt
    total_t > kstep && (Δt -= (total_t - kstep))

    # from there: kstep is replaced by this_step. LHE
    for i in 1:n
      if i == 1
        θ[i] += (inf - r_waterflow[i] - Ett[i]) * Δt / dz[i]
      else
        θ[i] += (r_waterflow[i-1] - r_waterflow[i] - Ett[i]) * Δt / dz[i]
      end
      θ[i] = clamp(θ[i], θ_vwp[i], θ_sat[i])
    end
  end

  for i in 1:n
    ice_ratio[i] *= θ_prev[i] / θ[i]
    ice_ratio[i] = min(1.0, ice_ratio[i])
  end
end


# Campbell 1974, Bonan 2019 Table 8.2
@fastmath function cal_ψ(θ::T, θ_sat::T, ψ_sat::T, b::T) where {T<:Real}
  ψ = ψ_sat * (θ / θ_sat)^(-b)
  max(ψ, ψ_sat)
end

@fastmath cal_K(θ::T, θ_sat::T, K_sat::T, b::T) where {T<:Real} =
  K_sat * (θ / θ_sat)^(2 * b + 3)

"""
[m s-1] -> 1000*[mm s-1] -> 1000*[kg m-2 s-1]
"""
# 如果流速过快，则减小时间步长
function guess_step(max_Fb)
  # this constraint is too large
  if max_Fb > 1.0e-5 # 864 mm/day
    dt = 1.0
  elseif max_Fb > 1.0e-6 # 86.4 mm/day
    dt = 30.0 # seconds
  else
    dt = 360.0
  end
  dt
end

# Function to calculate soil water uptake from a layer
"""
    Root Water Uptake

- `土壤蒸发`：仅发生在表层
- `植被蒸腾`：根据根系分布，耗水可能来自于土壤的每一层
"""
function Root_Water_Uptake(st::S, Trans_o::Float64, Trans_u::Float64, Evap_soil::Float64) where {
  S<:Union{StateBEPS,Soil}}

  Trans = Trans_o + Trans_u
  st.Ett[1] = Trans / ρ_w * st.dt[1] + Evap_soil / ρ_w
  for i in 2:st.n_layer
    st.Ett[i] = Trans / ρ_w * st.dt[i]
  end
end
