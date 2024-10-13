"""
o, u积雪的变化

- kstep: 360s, [s]
- snowrate: [m s-1]，注意`snowrate`的单位
"""
function snow_change(m_snow_pre::FT, snowrate::FT, kstep::FT,
  ρ_new_snow::FT, lai::FT, Ω::FT) where {FT<:Real}

  massMax_snow = 0.1 * lai
  areaMax_snow = 0.01 * lai

  τ = (1 - exp(-lai * Ω)) # 冠层截获的部分
  m_snow = m_snow_pre + snowrate * kstep * ρ_new_snow * τ

  perc_snow = m_snow / massMax_snow
  perc_snow = clamp(perc_snow, 0, 1)
  area_snow = perc_snow * areaMax_snow

  Δm_snow = m_snow - m_snow_pre
  return m_snow, perc_snow, area_snow, Δm_snow
end

"""
  snowpack_stage1_jl

# Arguments
- `m_snow`: m_snow
- `mw`: m_water
- `z_snow`: snow depth
- `z_water`: water depth

*reference variables*
- m_snow_pre
- m_snow
- perc_snow
- area_snow

# add an example of snowpack
"""
function snowpack_stage1_jl(Tair::Float64, prcp::Float64,
  lai_o::Float64, lai_u::Float64, Ω::Float64,
  m_snow_pre::Layer3{Float64},
  m_snow::Layer3{Float64},
  perc_snow::Layer3{Float64},
  area_snow::Layer2{Float64},
  z_snow::Float64,
  ρ_snow::Ref{Float64},
  albedo_v_snow::Ref{Float64}, albedo_n_snow::Ref{Float64})

  # m_snow_pre = Layer3(m_snow)
  massMax_snow_o = 0.1 * lai_o
  massMax_snow_u = 0.1 * lai_u

  # https://www.eoas.ubc.ca/courses/atsc113/snow/met_concepts/07-met_concepts/07b-newly-fallen-snow-density/
  ρ_new = 67.9 + 51.3 * exp(Tair / 2.6) # bug at here，新雪的密度
  # ρ_new = clamp(ρ_new, 50.0, 200.0) # 限制ρ_new的有效值域

  albedo_v_new = 0.94
  albedo_n_new = 0.8

  snowrate = Tair > 0 ? 0 : prcp * ρ_w / ρ_new

  mass2rate(Δm) = Δm / ρ_new / kstep # [kg m-2] -> [m s-1]
  cal_SnowPerc(m_snow, massMax_snow) = clamp(m_snow / massMax_snow, 0.0, 1.0)

  # kg2m 
  if Tair < 0
    snowrate_o = snowrate
    m_snow.o, perc_snow.o, area_snow.o, Δm_snow_o =
      snow_change(m_snow_pre.o, snowrate_o, kstep, ρ_new, lai_o, Ω)

    snowrate_u = max(0, snowrate_o - mass2rate(Δm_snow_o))
    m_snow.u, perc_snow.u, area_snow.u, Δm_snow_u =
      snow_change(m_snow_pre.u, snowrate_u, kstep, ρ_new, lai_u, Ω)

    snowrate_g = max(0.0, snowrate_u - mass2rate(Δm_snow_u))
    δ_zs = snowrate_g * kstep
  else
    snowrate_o = 0.0
    m_snow.o = m_snow_pre.o
    perc_snow.o = clamp(m_snow.o / massMax_snow_o, 0.0, 1.0)

    m_snow.u = m_snow_pre.u
    perc_snow.u = clamp(m_snow.u / massMax_snow_u, 0.0, 1.0)
    # area_snow.o = area_snow.o # area 不变
    # area_snow.u = area_snow.u
    δ_zs = 0.0
  end

  δ_zs = max(0.0, δ_zs)
  m_snow.g = max(0.0, m_snow_pre.g + δ_zs * ρ_new) # [kg m-2]

  if δ_zs > 0
    ρ_snow[] = (ρ_snow[] * z_snow + ρ_new * δ_zs) / (z_snow + δ_zs) # 计算混合密度
  else
    ρ_snow[] = (ρ_snow[] - 250) * exp(-0.001 * kstep / 3600.0) + 250.0
  end

  z_snow = m_snow.g > 0 ? m_snow.g / ρ_snow[] : 0.0
  perc_snow.g = min(m_snow.g / (0.05 * ρ_snow[]), 1.0) # [m]，认为雪深50cm时，perc_snow=1

  if snowrate_o > 0
    albedo_v_snow[] = (albedo_v_snow[] - 0.70) * exp(-0.005 * kstep / 3600) + 0.7
    albedo_n_snow[] = (albedo_n_snow[] - 0.42) * exp(-0.005 * kstep / 3600) + 0.42
  else
    albedo_v_snow[] = albedo_v_new
    albedo_n_snow[] = albedo_n_new
  end
  min(z_snow, 10.0) # 雪深过高，限制为10m即可
end


function snowpack_stage2_jl(evapo_snow_o::Float64, evapo_snow_u::Float64, m_snow::Layer3{Float64})
  # kstep::Float64 = kstep  # length of step
  m_snow.o = max(0.0, m_snow.o - evapo_snow_o * kstep)
  m_snow.u = max(0.0, m_snow.u - evapo_snow_u * kstep)
end


# 热量释放用于融雪（Tsnow -> 0）
# - pos: 融化
# - neg: 冻结
function cal_melt(z_snow::T, ρ_snow::T, Tsnow::T) where {T<:Real}
  dT = Tsnow - 0
  cp_ice = 2228.261         # J Kg-1 K-1
  λ_fusion = 3.34 * 1000000 # J Kg-1
  m = z_snow * ρ_snow   # [kg m-2]
  E = m * dT * cp_ice       # 当前的雪融化，需要这么多能量，
  return E / λ_fusion       # m_ice, -> kg
end

"""
Update Snow on the ground

It is assumed sublimation happens before the melting and freezing process.

> Note: 雪过深时，只有表层的温度(z=0.02)释放。这里存在bug。应该是
> `max(zs_sup*0.02, 0.02)`，而不是2%。

## Units
- mass: [kg m-2]
- depth: [m]
"""
function snowpack_stage3_jl(Tair::Float64, Tsnow::Float64, Tsnow_last::Float64, ρ_snow::Float64,
  z_snow::Float64, z_water::Float64, m_snow::Layer3{Float64})

  zs_sup = z_snow  # already considered sublimation
  ms_sup = m_snow.g

  Δm = cal_melt(zs_sup, ρ_snow, Tsnow) # [kg m-2]
  con_melt = Tsnow > 0 && Tsnow_last <= 0 && ms_sup > 0
  con_frozen = Tsnow <= 0 && Tsnow_last > 0 && z_water > 0

  ms_melt = 0.0
  mw_frozen = 0.0

  if zs_sup <= 0.02
    # case 1 depth of snow <0.02 m
    if Tair > 0 && zs_sup > 0
      ms_melt = min(Tair * 0.0075 * kstep / 3600 * 0.3, ms_sup)
    end
  elseif 0.02 < zs_sup <= 0.05
    # case 2 depth of snow > 0.02 < 0.05 m
    con_melt && (ms_melt = min(Δm, ms_sup))
    con_frozen && (mw_frozen = min(-Δm, z_water * ρ_w))
  elseif zs_sup > 0.05
    # _z = max(zs_sup*0.02, 0.02) # TODO: fix释放热量的深度
    # _Δm = cal_melt(_z, Tsnow)
    con_melt && (ms_melt = min(0.02Δm, ms_sup))
    con_frozen && (mw_frozen = min(-0.02Δm, z_water * ρ_w))
  end

  m_snow.g = max(0.0, m_snow.g - ms_melt + mw_frozen) # ground snow
  z_snow = max(0.0, zs_sup + (mw_frozen - ms_melt) / ρ_snow)
  z_water = max(0.0, z_water + (ms_melt - mw_frozen) / ρ_w)
  z_snow, z_water
end
