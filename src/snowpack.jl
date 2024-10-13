function snow_change(m_snow_pre::FT, snowrate::FT, kstep::FT,
  ρ_new_snow::FT, lai::FT, clumping::FT) where {FT<:Real}

  massMax_snow = 0.1 * lai
  areaMax_snow = 0.01 * lai

  m_snow = m_snow_pre + snowrate * kstep * ρ_new_snow * (1 - exp(-lai * clumping))
  perc_snow = m_snow / massMax_snow
  perc_snow = clamp(perc_snow, 0, 1)
  area_snow = perc_snow * areaMax_snow
  massStep_snow = m_snow - m_snow_pre

  return m_snow, perc_snow, area_snow, massStep_snow
end

"""
  snowpack_stage1_jl

# Arguments

- `m_snow`: m_snow
- `mw`: m_water
- `depth_snow`: snow depth
- `depth_water`: water depth

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
  depth_snow::Float64,
  ρ_snow::Ref{Float64},
  albedo_v_snow::Ref{Float64}, albedo_n_snow::Ref{Float64})

  # m_snow_pre = Layer3(m_snow)
  massMax_snow_o = 0.1 * lai_o
  massMax_snow_u = 0.1 * lai_u

  ρ_new_snow = 67.9 + 51.3 * exp(Tair / 2.6) # bug at here
  albedo_v_Newsnow = 0.94
  albedo_n_Newsnow = 0.8

  snowrate = Tair > 0 ? 0 : prcp * ρ_w / ρ_new_snow
  snowrate_o = 0.0
  
  if Tair < 0
    snowrate_o = snowrate
    m_snow.o, perc_snow.o, area_snow.o, massStep_snow_o =
      snow_change(m_snow_pre.o, snowrate_o, kstep, ρ_new_snow, lai_o, Ω)

    snowrate_u = max(0, snowrate_o - (massStep_snow_o) / ρ_new_snow / kstep)
    m_snow.u, perc_snow.u, area_snow.u, massStep_snow_u =
      snow_change(m_snow_pre.u, snowrate_u, kstep, ρ_new_snow, lai_u, Ω)

    snowrate_g = max(0.0, snowrate_u - (massStep_snow_u) / ρ_new_snow / kstep)
    δ_zs = snowrate_g * kstep
  else
    m_snow.o = m_snow_pre.o
    perc_snow.o = clamp(m_snow.o / massMax_snow_o, 0.0, 1.0)

    m_snow.u = m_snow_pre.u
    perc_snow.u = clamp(m_snow.u / massMax_snow_u, 0.0, 1.0)
    # area_snow.o = area_snow.o # area 不变
    # area_snow.u = area_snow.u
    δ_zs = 0.0
  end

  δ_zs = max(0.0, δ_zs)
  m_snow.g = max(0.0, m_snow_pre.g + δ_zs * ρ_new_snow)

  if δ_zs > 0
    ρ_snow[] = (ρ_snow[] * depth_snow + ρ_new_snow * δ_zs) / (depth_snow + δ_zs)
  else
    ρ_snow[] = (ρ_snow[] - 250) * exp(-0.001 * kstep / 3600.0) + 250.0
  end

  depth_snow = m_snow.g > 0 ? m_snow.g / ρ_snow[] : 0.0
  perc_snow.g = min(m_snow.g / (0.05 * ρ_snow[]), 1.0)

  if snowrate_o > 0
    albedo_v_snow[] = (albedo_v_snow[] - 0.70) * exp(-0.005 * kstep / 3600) + 0.7
    albedo_n_snow[] = (albedo_n_snow[] - 0.42) * exp(-0.005 * kstep / 3600) + 0.42
  else
    albedo_v_snow[] = albedo_v_Newsnow
    albedo_n_snow[] = albedo_n_Newsnow
  end

  min(depth_snow, 100.0)
end


function snowpack_stage2_jl(evapo_snow_o::Float64, evapo_snow_u::Float64, m_snow::Layer3{Float64})
  # kstep::Float64 = kstep  # length of step
  m_snow.o = max(0.0, m_snow.o - evapo_snow_o * kstep)
  m_snow.u = max(0.0, m_snow.u - evapo_snow_u * kstep)
end


function snowpack_stage3_jl(Tair::Float64, Tsnow::Float64, Tsnow_last::Float64, ρ_snow::Float64,
  depth_snow::Float64, depth_water::Float64, m_snow::Layer3{Float64})

  # kstep = kstep  # kstep in model
  # it is assumed sublimation happens before the melting and freezing process
  zs_sup = depth_snow  # this snow depth has already considered sublimation
  ms_sup = m_snow.g

  # parameters
  cp_ice = 2228.261
  latent_fusion = 3.34 * 1000000 # J Kg-1
  ρ_w = 1025 # Kg m-3

  # simulate snow melt and freeze process
  ms_melt = 0.0
  mw_frozen = 0.0

  # case 1 depth of snow <0.02 m
  if zs_sup <= 0.02
    if Tair > 0 && zs_sup > 0
      ms_melt = Tair * 0.0075 * kstep / 3600 * 0.3
      ms_melt = min(ms_sup, ms_melt)  # the amount of melted snow could not be larger than supply
    else
      ms_melt = 0.0
    end
    # case 2 depth of snow > 0.02 < 0.05 m
  elseif 0.02 < zs_sup <= 0.05
    if Tsnow > 0 && Tsnow_last < 0 && ms_sup > 0  # melting
      ms_melt = Tsnow * cp_ice * ρ_snow * zs_sup / latent_fusion
      ms_melt = min(ms_sup, ms_melt)  # the amount of melted snow could not be larger than supply
    else
      ms_melt = 0.0
    end

    if Tsnow <= 0 && Tsnow_last > 0 && depth_water > 0  # freezing
      mw_frozen = (0 - Tsnow) * cp_ice * ρ_snow * zs_sup / latent_fusion
      mw_frozen = max(mw_frozen, depth_water * ρ_w)
    else
      mw_frozen = 0.0
    end
    # case 3 depth of snow > 0.05 m
  elseif zs_sup > 0.05
    if Tsnow > 0 && Tsnow_last <= 0 && ms_sup > 0  # melting
      ms_melt = Tsnow * cp_ice * ρ_snow * zs_sup * 0.02 / latent_fusion
      ms_melt = min(ms_sup, ms_melt)  # the amount of melted snow could not be larger than supply
    else
      ms_melt = 0.0
    end

    if Tsnow <= 0 && Tsnow_last > 0 && depth_water > 0  # freezing
      mw_frozen = (0 - Tsnow) * cp_ice * ρ_snow * zs_sup * 0.02 / latent_fusion
      mw_frozen = max(mw_frozen, depth_water * ρ_w)
    else
      mw_frozen = 0.0
    end
  end

  # δ in mass of snow on ground
  m_snow.g = max(0.0, m_snow.g - ms_melt + mw_frozen)

  # δ of depth in snow
  melt_zs = ms_melt / ρ_snow
  frozen_zs = mw_frozen / ρ_snow
  depth_snow = max(0.0, zs_sup - melt_zs + frozen_zs)

  # δ of depth in water
  melt_zw = ms_melt / ρ_w
  frozen_zw = mw_frozen / ρ_w
  depth_water = max(0.0, depth_water + melt_zw - frozen_zw)

  depth_snow, depth_water
end
