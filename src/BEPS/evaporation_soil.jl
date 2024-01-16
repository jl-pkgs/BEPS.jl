function evaporation_soil_jl(Tair::FT, Tg::FT, RH::FT, Rn_g::FT, Gheat_g::FT,
  percent_snow_g::Ref{FT}, depth_water::Ref{FT}, depth_snow::Ref{FT}, mass_water_g::Ref{FT}, mass_snow_g::Ref{FT}, # Ref{FT}
  density_snow::FT, swc_g::FT, porosity_g::FT) where {FT<:Real}
  # evapo_soil::Ref{FT}, evapo_water_g::Ref{FT}, evapo_snow_g::Ref{FT}

  met = meteo_pack_jl(Tg, RH)
  ρ_a = met.rho_a
  cp = met.cp
  vpd = met.VPD
  Δ = met.slope
  γ = met.gamma

  λ_water = (2.501 - 0.00237 * Tair) * 1000000
  λ_snow = 2.83 * 1000000
  ρ_water = 1025.0

  Gwater_g = 1.0 / (4.0 * exp(8.2 - 4.2 * swc_g / porosity_g))
  length_step = kstep

  if depth_snow[] > 0.02
    percent_snow_g[] = 1.0
  else
    percent_snow_g[] = mass_snow_g[] / (0.025 * density_snow)
  end

  percent_snow_g[] = max(percent_snow_g[], 0.0)
  percent_snow_g[] = min(percent_snow_g[], 1.0)

  
  if depth_water[] > 0.0 && depth_snow[] == 0.0
    evapo_water_g = 1.0 / (λ_water) * (Δ * (Rn_g * 0.8 - 0) + ρ_a * cp * vpd * Gheat_g) / (Δ + γ * (1 + (Gheat_g) / 0.01))
  else
    evapo_water_g = 0.0
  end

  evapo_water_g = max(-0.002 / length_step, evapo_water_g)
  if evapo_water_g > 0
    evapo_water_g = min(evapo_water_g, (depth_water[] * ρ_water) / length_step)
  end

  depth_water[] = depth_water[] - (evapo_water_g / ρ_water) * length_step
  depth_water[] = max(0, depth_water[])
  mass_water_g[] = mass_water_g[] - evapo_water_g * length_step

  if depth_snow[] > 0.0
    evapo_snow_g = 1 / (λ_snow) * (Δ * (Rn_g * 0.8 - 0) + ρ_a * cp * vpd * Gheat_g) / (Δ + γ * (1 + (Gheat_g) / 0.01)) * (percent_snow_g[])
  else
    evapo_snow_g = 0.0
  end

  evapo_snow_g = max(-0.002 / length_step, evapo_snow_g)
  if evapo_snow_g > 0
    evapo_snow_g = min(evapo_snow_g, mass_snow_g[] / length_step)
  end

  mass_snow_g[] = mass_snow_g[] - evapo_snow_g * length_step
  mass_snow_g[] = max(mass_snow_g[], 0)

  if mass_snow_g[] > 0.0
    depth_snow[] = depth_snow[] - (evapo_snow_g / density_snow) * length_step
  else
    depth_snow[] = 0.0
  end

  if depth_water[] > 0.0 || depth_snow[] > 0.0
    evapo_soil = 0.0
  else
    evapo_soil = (1.0 - (percent_snow_g[])) * 1 / (λ_water) * (Δ * (Rn_g - 0) + ρ_a * cp * vpd * Gheat_g) / (Δ + γ * (1 + Gheat_g / Gwater_g))
    evapo_soil = max(0.0, evapo_soil)
  end

  evapo_soil, evapo_water_g, evapo_snow_g
  # return results
end
