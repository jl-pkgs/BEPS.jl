abstract type AbstractAirLayer{FT} end


"""
    mutable struct AirLayer{FT}

Struct to store environmental conditions in each air layer corresponds to one canopy layer.

# Fields
$(TYPEDFIELDS)

> Copied from Land.jl
"""
@with_kw mutable struct AirLayer{FT<:AbstractFloat} <: AbstractAirLayer{FT}
  "Air temperature `[K]`"
  Tair::FT
  
  "Air density `[kg m⁻³]`"
  rho_a::FT
  
  "Specific heat of air `[J kg⁻¹ K⁻¹]`"
  Cp_ca::FT

  "Vapor pressure deficit `[Pa]`"
  VPD::FT

  "Psychrometric constant `[Pa K⁻¹]`"
  γ::FT

  "Slope of Saturation vapor pressure es `[Pa K⁻¹]`"
  Δ::FT

  "Retive humility `[%]`"
  RH::FT

  "Wind speed `[m s⁻¹]`"
  wind::FT = FT(2) 
end

# args = (Tair, slope, γ, VPD_air, Cp_ca)

export AirLayer
