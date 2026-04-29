"""
Inverse diagnostic tools for aerodynamic resistance and canopy conductance.

Provides three complementary approaches:
- Method A: direct ra inversion from observed SH + LST
- Method B: Penman-Monteith inversion of Gc from observed LE
- Helper:   stability class labeling for stratified analysis
"""


"""
    ra_from_flux(SH_obs, Tc, Ta; ρₐ=1.292, cp=1010.0) -> ra [s/m]

Method A: directly invert aerodynamic resistance from observed sensible heat
flux and canopy (or surface) temperature.

# Physics
```
SH = ρₐ cp (Tc - Ta) / ra  →  ra = ρₐ cp (Tc - Ta) / SH
```

# Arguments
- `SH_obs` : observed sensible heat flux [W/m²]
- `Tc`     : canopy / surface radiometric temperature [°C]
- `Ta`     : air temperature at reference height [°C]
- `ρₐ`     : air density [kg/m³]  (default 1.292)
- `cp`     : specific heat of air [J/kg/K]  (default 1010.0)

# Returns
`ra` in [s/m]; returns `Inf` when `|SH_obs| < 1 W/m²` (near-zero flux
condition where the inversion is physically unreliable — such points should be
filtered before further analysis).
"""
function ra_from_flux(SH_obs::T, Tc::T, Ta::T;
  ρₐ::T=T(1.292), cp::T=T(1010.0))::T where {T<:Real}
  abs(SH_obs) < 1.0 && return T(Inf)
  ρₐ * cp * (Tc - Ta) / SH_obs
end


"""
    Gc_penman(LE_obs, Ga, Rn, G_soil, Ta, RH) -> Gc [m/s]

Method B: invert bulk canopy conductance from observed latent heat flux using
the Penman-Monteith equation rearranged for Gc.

# Physics
PM equation solved for Gc:
```
LE = [Δ(Rn-G) + ρₐ cp VPD Ga] / [Δ + γ(1 + Ga/Gc)]
  →  Gc = γ LE Ga / [Δ(Rn-G) + ρₐ cp VPD Ga - (Δ+γ) LE]
```

# Arguments
- `LE_obs` : observed latent heat flux [W/m²]
- `Ga`     : aerodynamic conductance for heat [m/s]
- `Rn`     : net radiation [W/m²]
- `G_soil` : soil heat flux [W/m²]
- `Ta`     : air temperature [°C]
- `RH`     : relative humidity [%]

# Returns
`Gc` in [m/s]; returns `NaN` when denominator < 1e-3 (near-saturated or
energy-limited regime where PM inversion is ill-conditioned — filter `isnan`
results before statistical analysis).
"""
function Gc_penman(LE_obs::T, Ga::T, Rn::T, G_soil::T,
  Ta::T, RH::T)::T where {T<:Real}
  met = meteo_pack_jl(Ta, RH)
  (; Δ, γ, ρₐ, cp, VPD) = met

  denom = Δ * (Rn - G_soil) + ρₐ * cp * VPD * Ga - (Δ + γ) * LE_obs
  abs(denom) < 1e-3 && return T(NaN)
  γ * LE_obs * Ga / denom
end


"""
    stability_class(ξ) -> Symbol

Classify atmospheric stability from the Monin-Obukhov stability parameter
ξ = (z-d)/L.

| ξ          | Class          |
|:-----------|:---------------|
| < -0.5     | :unstable      |
| -0.5 to 0.5| :neutral       |
| > 0.5      | :stable        |
"""
function stability_class(ξ::Real)::Symbol
  ξ < -0.5 && return :unstable
  ξ > 0.5  && return :stable
  return :neutral
end
