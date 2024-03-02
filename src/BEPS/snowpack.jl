"""
  snowpack_stage1_jl

- `ms`: mass_snow
- `mw`: mass_water
- `d_s`: snow depth
- `d_w`: water depth
"""
function snowpack_stage1_jl(temp_air::Float64, prcp::Float64,
  # ms_o_last::Float64, ms_u_last::Float64, ms_g_last::Float64, 
  ms::CanopyLayer{Float64},
  # ms_o::Ref{Float64}, ms_u::Ref{Float64}, ms_g::Ref{Float64},
  lai_o::Float64, lai_u::Float64, clumping::Float64, area_snow_o::Ref{Float64}, area_snow_u::Ref{Float64},
  # perc_snow_o::Ref{Float64}, perc_snow_u::Ref{Float64}, perc_snow_g::Ref{Float64},
  ρ_snow::Ref{Float64}, d_s, albedo_v_snow::Ref{Float64}, albedo_n_snow::Ref{Float64})

  massMax_snow_o = 0.1 * lai_o
  massMax_snow_u = 0.1 * lai_u
  areaMax_snow_o = lai_o * 0.01
  areaMax_snow_u = lai_u * 0.01

  ρ_new_snow = 67.9 + 51.3 * exp(temp_air / 2.6)
  albedo_v_Newsnow = 0.94
  albedo_n_Newsnow = 0.8

  ms_o_last = ms.o
  ms_u_last = ms.u
  ms_g_last = ms.g

  snowrate = temp_air > 0 ? 0 : prcp * ρ_w / ρ_new_snow
  snowrate_o = 0.0

  if temp_air < 0
    snowrate_o = snowrate
    ms.o = ms_o_last + snowrate_o * kstep * ρ_new_snow * (1 - exp(-lai_o * clumping))
    perc_snow_o = ms.o / massMax_snow_o
    perc_snow_o = clamp(perc_snow_o, 0, 1)
    area_snow_o[] = perc_snow_o * areaMax_snow_o
    massStep_snow_o = ms.o - ms_o_last

    snowrate_u = max(0, snowrate_o - (massStep_snow_o) / ρ_new_snow / kstep)
    ms.u = ms_u_last + snowrate_u * kstep * ρ_new_snow * (1 - exp(-lai_u * clumping))
    perc_snow_u = ms.u / massMax_snow_u
    perc_snow_u = clamp(perc_snow_u, 0, 1)
    area_snow_u[] = perc_snow_u * areaMax_snow_u
    massStep_snow_u = ms.u - ms_u_last

    snowrate_g = max(0, snowrate_u - (massStep_snow_u) / ρ_new_snow / kstep)
    δ_d_s = snowrate_g * kstep
  else
    ms.o = ms_o_last
    perc_snow_o = clamp(ms.o / massMax_snow_o, 0, 1)
    area_snow_o[] = area_snow_o[]

    ms.u = ms_u_last
    perc_snow_u = ms.u / massMax_snow_u
    perc_snow_u = clamp(perc_snow_u, 0, 1)
    area_snow_u[] = area_snow_u[]

    δ_d_s = 0
  end

  δ_d_s = max(0, δ_d_s)
  ms.g = max(0, ms_g_last + δ_d_s * ρ_new_snow)

  if δ_d_s > 0
    ρ_snow[] = (ρ_snow[] * d_s + ρ_new_snow * δ_d_s) / (d_s + δ_d_s)
  else
    ρ_snow[] = (ρ_snow[] - 250) * exp(-0.001 * kstep / 3600.0) + 250.0
  end

  d_s = ms.g > 0 ? ms.g / ρ_snow[] : 0.0
  perc_snow_g = min(ms.g / (0.05 * ρ_snow[]), 1.0)

  if snowrate_o > 0
    albedo_v_snow[] = (albedo_v_snow[] - 0.70) * exp(-0.005 * kstep / 3600) + 0.7
    albedo_n_snow[] = (albedo_n_snow[] - 0.42) * exp(-0.005 * kstep / 3600) + 0.42
  else
    albedo_v_snow[] = albedo_v_Newsnow
    albedo_n_snow[] = albedo_n_Newsnow
  end

  perc_snow = CanopyLayer(perc_snow_o, perc_snow_u, perc_snow_g)
  perc_snow, d_s
end


function snowpack_stage2_jl(evapo_snow_o::Float64, evapo_snow_u::Float64, ms::CanopyLayer{Float64})
  # ms_o::Ref{Float64}, ms_u::Ref{Float64})

  # kstep::Float64 = kstep  # length of step
  ms.o = max(0, ms.o - evapo_snow_o * kstep)
  ms.u = max(0, ms.u - evapo_snow_u * kstep)
end


function snowpack_stage3_jl(temp_air::Float64, temp_snow::Float64, temp_snow_last::Float64, ρ_snow::Float64,
  d_s, d_w,
  # d_s::Ref{Float64}, d_w::Ref{Float64}, 
  ms::CanopyLayer{Float64})
  # ms_g::Ref{Float64})

  # kstep = kstep  # kstep in model
  # it is assumed sublimation happens before the melting and freezing process
  d_s_sup = d_s  # this snow depth has already considered sublimation
  ms_sup = ms.g

  # parameters
  cp_ice = 2228.261
  latent_fusion = 3.34 * 1000000 # J Kg-1
  ρ_w = 1025 # Kg m-3

  # simulate snow melt and freeze process
  ms_melt = 0.0
  mw_frozen = 0.0

  # case 1 depth of snow <0.02 m
  if d_s_sup <= 0.02
    if temp_air > 0 && d_s_sup > 0
      ms_melt = temp_air * 0.0075 * kstep / 3600 * 0.3
      ms_melt = min(ms_sup, ms_melt)  # the amount of melted snow could not be larger than supply
    else
      ms_melt = 0
    end
    # case 2 depth of snow > 0.02 < 0.05 m
  elseif d_s_sup > 0.02 && d_s_sup <= 0.05
    if temp_snow > 0 && temp_snow_last < 0 && ms_sup > 0  # melting
      ms_melt = temp_snow * cp_ice * ρ_snow * d_s_sup / latent_fusion
      ms_melt = min(ms_sup, ms_melt)  # the amount of melted snow could not be larger than supply
    else
      ms_melt = 0
    end

    if temp_snow <= 0 && temp_snow_last > 0 && d_w > 0  # freezing
      mw_frozen = (0 - temp_snow) * cp_ice * ρ_snow * d_s_sup / latent_fusion
      mw_frozen = max(mw_frozen, d_w * ρ_w)
    else
      mw_frozen = 0
    end
    # case 3 depth of snow > 0.05 m
  elseif d_s_sup > 0.05
    if temp_snow > 0 && temp_snow_last <= 0 && ms_sup > 0  # melting
      ms_melt = temp_snow * cp_ice * ρ_snow * d_s_sup * 0.02 / latent_fusion
      ms_melt = min(ms_sup, ms_melt)  # the amount of melted snow could not be larger than supply
    else
      ms_melt = 0
    end

    if temp_snow <= 0 && temp_snow_last > 0 && d_w > 0  # freezing
      mw_frozen = (0 - temp_snow) * cp_ice * ρ_snow * d_s_sup * 0.02 / latent_fusion
      mw_frozen = max(mw_frozen, d_w * ρ_w)
    else
      mw_frozen = 0
    end
  end

  # δ in mass of snow on ground
  ms.g = ms.g - ms_melt + mw_frozen
  ms.g = max(0, ms.g)

  # δ of depth in snow
  melt_d_s = ms_melt / ρ_snow
  frozen_d_s = mw_frozen / ρ_snow
  d_s = d_s_sup - melt_d_s + frozen_d_s
  d_s = max(0, d_s)

  # δ of depth in water
  melt_d_w = ms_melt / ρ_w
  frozen_d_w = mw_frozen / ρ_w
  d_w = d_w + melt_d_w - frozen_d_w
  d_w = max(0, d_w)

  d_s, d_w
end
