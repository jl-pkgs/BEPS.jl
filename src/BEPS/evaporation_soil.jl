"""

# Return
Ew: evaporation from water
Es: evaporation from soil
"""
function evaporation_soil_jl(Tair::FT, Tg::FT, RH::FT, Rn_g::FT, Gheat_g::FT,
  perc_snow_g::Ref{FT}, depth_water::Ref{FT}, depth_snow::Ref{FT}, mass_water_g::Ref{FT}, mass_snow_g::Ref{FT},
  ρ_snow::FT, swc_g::FT, porosity_g::FT) where {FT<:Real}
  
  met = meteo_pack_jl(Tg, RH)
  ρ_a = met.ρ_a
  cp = met.cp
  vpd = met.VPD
  Δ = met.slope
  γ = met.gamma

  λ_water = (2.501 - 0.00237 * Tair) * 1000000
  
  Gwater_g = 1.0 / (4.0 * exp(8.2 - 4.2 * swc_g / porosity_g))
  length_step = kstep

  perc_snow_g[] = depth_snow[] > 0.02 ? 1.0 : mass_snow_g[] / (0.025 * ρ_snow)
  perc_snow_g[] = clamp(perc_snow_g[], 0.0, 1.0)

  if depth_water[] > 0.0 && depth_snow[] == 0.0
    Ew_g = 1.0 / (λ_water) * (Δ * (Rn_g * 0.8 - 0) + ρ_a * cp * vpd * Gheat_g) / (Δ + γ * (1 + (Gheat_g) / 0.01))
  else
    Ew_g = 0.0
  end

  Ew_g = max(-0.002 / length_step, Ew_g)
  if Ew_g > 0.0
    Ew_g = min(Ew_g, (depth_water[] * ρ_w) / length_step)
  end

  depth_water[] = depth_water[] - (Ew_g / ρ_w) * length_step
  depth_water[] = max(0, depth_water[])
  mass_water_g[] = mass_water_g[] - Ew_g * length_step

  if depth_snow[] > 0.0
    Esoil_g = 1 / (λ_snow) * (Δ * (Rn_g * 0.8 - 0) + ρ_a * cp * vpd * Gheat_g) / (Δ + γ * (1 + (Gheat_g) / 0.01)) * perc_snow_g[]
  else
    Esoil_g = 0.0
  end

  Esoil_g = max(-0.002 / length_step, Esoil_g)
  if Esoil_g > 0.0
    Esoil_g = min(Esoil_g, mass_snow_g[] / length_step)
  end

  mass_snow_g[] = max(mass_snow_g[] - Esoil_g * length_step, 0)
  depth_snow[] = mass_snow_g[] > 0.0 ? depth_snow[] - (Esoil_g / ρ_snow) * length_step : 0.0

  if depth_water[] > 0.0 || depth_snow[] > 0.0
    Es = 0.0
  else
    Es = (1.0 - (perc_snow_g[])) * 1 / (λ_water) * (Δ * (Rn_g - 0) + ρ_a * cp * vpd * Gheat_g) / (Δ + γ * (1 + Gheat_g / Gwater_g))
    Es = max(0.0, Es)
  end

  Es, Ew_g, Esoil_g
end
