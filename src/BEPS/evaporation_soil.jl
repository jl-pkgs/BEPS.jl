"""

# Return
Ewater: evaporation from water
Esoil: evaporation from soil
"""
function evaporation_soil_jl(Tair::FT, Tg::FT, RH::FT, Rn_g::FT, Gheat_g::FT,
  perc_snow_g::Ref{FT}, depth_water::Ref{FT}, depth_snow::Ref{FT},
  mass_water_g::Ref{FT}, mass_snow_g::Ref{FT},
  ρ_snow::FT, swc_g::FT, porosity_g::FT) where {FT<:Real}

  met = meteo_pack_jl(Tg, RH)
  (; ρₐ, cp, VPD, Δ, γ) = met
  λ = cal_lambda(Tair) #

  Gwater_g = 1.0 / (4.0 * exp(8.2 - 4.2 * swc_g / porosity_g))
  # kstep = kstep # 360s, 6min

  perc_snow_g[] = depth_snow[] > 0.02 ? 1.0 : mass_snow_g[] / (0.025 * ρ_snow)
  perc_snow_g[] = clamp(perc_snow_g[], 0.0, 1.0)

  if depth_water[] > 0.0 && depth_snow[] == 0.0
    Ewater_g = 1.0 / λ * (Δ * (Rn_g * 0.8 - 0) + ρₐ * cp * VPD * Gheat_g) /
               (Δ + γ * (1 + (Gheat_g) / 0.01))
  else
    Ewater_g = 0.0
  end

  Ewater_g = max(-0.002 / kstep, Ewater_g)
  if Ewater_g > 0.0
    Ewater_g = min(Ewater_g, (depth_water[] * ρ_w) / kstep)
  end

  depth_water[] = depth_water[] - (Ewater_g / ρ_w) * kstep
  depth_water[] = max(0, depth_water[])
  mass_water_g[] = mass_water_g[] - Ewater_g * kstep

  if depth_snow[] > 0.0
    Esoil_g = 1 / λ_snow * (Δ * (Rn_g * 0.8 - 0) + ρₐ * cp * VPD * Gheat_g) /
              (Δ + γ * (1 + (Gheat_g) / 0.01)) * perc_snow_g[]
  else
    Esoil_g = 0.0
  end

  Esoil_g = max(-0.002 / kstep, Esoil_g)
  if Esoil_g > 0.0
    Esoil_g = min(Esoil_g, mass_snow_g[] / kstep)
  end

  mass_snow_g[] = max(mass_snow_g[] - Esoil_g * kstep, 0)
  depth_snow[] = mass_snow_g[] > 0.0 ? depth_snow[] - (Esoil_g / ρ_snow) * kstep : 0.0

  if depth_water[] > 0.0 || depth_snow[] > 0.0
    Esoil = 0.0
  else
    Esoil = (1.0 - (perc_snow_g[])) * 1 / (λ) * (Δ * (Rn_g - 0) + ρₐ * cp * VPD * Gheat_g) /
            (Δ + γ * (1 + Gheat_g / Gwater_g))
    Esoil = max(0.0, Esoil)
  end

  Esoil, Ewater_g, Esoil_g
end
