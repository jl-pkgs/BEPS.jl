function snowpack_stage1_jl(temp_air::Float64, prcp::Float64,
  # mass_snow_o_last::Float64, mass_snow_u_last::Float64, mass_snow_g_last::Float64, 
  mass_snow_o::Ref{Float64}, mass_snow_u::Ref{Float64}, mass_snow_g::Ref{Float64},
  lai_o::Float64, lai_u::Float64, clumping::Float64, area_snow_o::Ref{Float64}, area_snow_u::Ref{Float64},
  # perc_snow_o::Ref{Float64}, perc_snow_u::Ref{Float64}, perc_snow_g::Ref{Float64},
  ρ_snow::Ref{Float64}, depth_snow::Ref{Float64}, albedo_v_snow::Ref{Float64}, albedo_n_snow::Ref{Float64})

  massMax_snow_o = 0.1 * lai_o
  massMax_snow_u = 0.1 * lai_u
  areaMax_snow_o = lai_o * 0.01
  areaMax_snow_u = lai_u * 0.01

  length_step = kstep

  ρ_new_snow = 67.9 + 51.3 * exp(temp_air / 2.6)
  ρ_water = 1025
  albedo_v_Newsnow = 0.94
  albedo_n_Newsnow = 0.8

  mass_snow_o_last = mass_snow_o[]
  mass_snow_u_last = mass_snow_u[]
  mass_snow_g_last = mass_snow_g[]

  snowrate = temp_air > 0 ? 0 : prcp * ρ_water / ρ_new_snow
  snowrate_o = 0.0

  if temp_air < 0
    snowrate_o = snowrate
    mass_snow_o[] = mass_snow_o_last + snowrate_o * length_step * ρ_new_snow * (1 - exp(-lai_o * clumping))
    perc_snow_o = mass_snow_o[] / massMax_snow_o
    perc_snow_o = clamp(perc_snow_o, 0, 1)
    area_snow_o[] = perc_snow_o * areaMax_snow_o
    massStep_snow_o = mass_snow_o[] - mass_snow_o_last

    snowrate_u = max(0, snowrate_o - (massStep_snow_o) / ρ_new_snow / length_step)
    mass_snow_u[] = mass_snow_u_last + snowrate_u * length_step * ρ_new_snow * (1 - exp(-lai_u * clumping))
    perc_snow_u = mass_snow_u[] / massMax_snow_u
    perc_snow_u = clamp(perc_snow_u, 0, 1)
    area_snow_u[] = perc_snow_u * areaMax_snow_u
    massStep_snow_u = mass_snow_u[] - mass_snow_u_last

    snowrate_g = max(0, snowrate_u - (massStep_snow_u) / ρ_new_snow / length_step)
    change_depth_snow = snowrate_g * length_step
  else
    mass_snow_o[] = mass_snow_o_last
    perc_snow_o = clamp(mass_snow_o[] / massMax_snow_o, 0, 1)
    area_snow_o[] = area_snow_o[]

    mass_snow_u[] = mass_snow_u_last
    perc_snow_u = mass_snow_u[] / massMax_snow_u
    perc_snow_u = clamp(perc_snow_u, 0, 1)
    area_snow_u[] = area_snow_u[]

    change_depth_snow = 0
  end

  change_depth_snow = max(0, change_depth_snow)
  mass_snow_g[] = mass_snow_g_last + change_depth_snow * ρ_new_snow
  mass_snow_g[] = max(0, mass_snow_g[])

  if change_depth_snow > 0
    ρ_snow[] = (ρ_snow[] * depth_snow[] + ρ_new_snow * change_depth_snow) / (depth_snow[] + change_depth_snow)
  else
    ρ_snow[] = (ρ_snow[] - 250) * exp(-0.001 * length_step / 3600.0) + 250.0
  end

  depth_snow[] = mass_snow_g[] > 0 ? mass_snow_g[] / ρ_snow[] : 0.0
  perc_snow_g = min(mass_snow_g[] / (0.05 * ρ_snow[]), 1.0)

  if snowrate_o > 0
    albedo_v_snow[] = (albedo_v_snow[] - 0.70) * exp(-0.005 * length_step / 3600) + 0.7
    albedo_n_snow[] = (albedo_n_snow[] - 0.42) * exp(-0.005 * length_step / 3600) + 0.42
  else
    albedo_v_snow[] = albedo_v_Newsnow
    albedo_n_snow[] = albedo_n_Newsnow
  end

  perc_snow_o, perc_snow_u, perc_snow_g
end


function snowpack_stage2_jl(evapo_snow_o::Float64, evapo_snow_u::Float64,
  mass_snow_o::Ref{Float64}, mass_snow_u::Ref{Float64})

  length_step::Float64 = kstep  # length of step

  mass_snow_o[] = max(0, mass_snow_o[] - evapo_snow_o * length_step)
  mass_snow_u[] = max(0, mass_snow_u[] - evapo_snow_u * length_step)
end


function snowpack_stage3_jl(temp_air::Float64, temp_snow::Float64, temp_snow_last::Float64, ρ_snow::Float64,
  depth_snow::Ref{Float64}, depth_water::Ref{Float64}, mass_snow_g::Ref{Float64})

  length_step = kstep  # length_step in model

  # it is assumed sublimation happens before the melting and freezing process
  depth_snow_sup = depth_snow[]  # this snow depth has already considered sublimation
  mass_snow_sup = mass_snow_g[]

  # parameters
  cp_ice = 2228.261
  latent_fusion = 3.34 * 1000000 # J Kg-1
  ρ_water = 1025 # Kg m-3

  # simulate snow melt and freeze process
  mass_snow_melt = 0.0
  mass_water_frozen = 0.0

  # case 1 depth of snow <0.02 m
  if depth_snow_sup <= 0.02
    if temp_air > 0 && depth_snow_sup > 0
      mass_snow_melt = temp_air * 0.0075 * length_step / 3600 * 0.3
      mass_snow_melt = min(mass_snow_sup, mass_snow_melt)  # the amount of melted snow could not be larger than supply
    else
      mass_snow_melt = 0
    end
    # case 2 depth of snow > 0.02 < 0.05 m
  elseif depth_snow_sup > 0.02 && depth_snow_sup <= 0.05
    if temp_snow > 0 && temp_snow_last < 0 && mass_snow_sup > 0  # melting
      mass_snow_melt = temp_snow * cp_ice * ρ_snow * depth_snow_sup / latent_fusion
      mass_snow_melt = min(mass_snow_sup, mass_snow_melt)  # the amount of melted snow could not be larger than supply
    else
      mass_snow_melt = 0
    end

    if temp_snow <= 0 && temp_snow_last > 0 && depth_water[] > 0  # freezing
      mass_water_frozen = (0 - temp_snow) * cp_ice * ρ_snow * depth_snow_sup / latent_fusion
      mass_water_frozen = max(mass_water_frozen, depth_water[] * ρ_water)
    else
      mass_water_frozen = 0
    end
    # case 3 depth of snow > 0.05 m
  elseif depth_snow_sup > 0.05
    if temp_snow > 0 && temp_snow_last <= 0 && mass_snow_sup > 0  # melting
      mass_snow_melt = temp_snow * cp_ice * ρ_snow * depth_snow_sup * 0.02 / latent_fusion
      mass_snow_melt = min(mass_snow_sup, mass_snow_melt)  # the amount of melted snow could not be larger than supply
    else
      mass_snow_melt = 0
    end

    if temp_snow <= 0 && temp_snow_last > 0 && depth_water[] > 0  # freezing
      mass_water_frozen = (0 - temp_snow) * cp_ice * ρ_snow * depth_snow_sup * 0.02 / latent_fusion
      mass_water_frozen = max(mass_water_frozen, depth_water[] * ρ_water)
    else
      mass_water_frozen = 0
    end
  end

  # change in mass of snow on ground
  mass_snow_g[] = mass_snow_g[] - mass_snow_melt + mass_water_frozen
  mass_snow_g[] = max(0, mass_snow_g[])

  # change of depth in snow
  melt_depth_snow = mass_snow_melt / ρ_snow
  frozen_depth_snow = mass_water_frozen / ρ_snow
  depth_snow[] = depth_snow_sup - melt_depth_snow + frozen_depth_snow
  depth_snow[] = max(0, depth_snow[])

  # change of depth in water
  melt_depth_water = mass_snow_melt / ρ_water
  frozen_depth_water = mass_water_frozen / ρ_water
  depth_water[] = depth_water[] + melt_depth_water - frozen_depth_water
  depth_water[] = max(0, depth_water[])
end
