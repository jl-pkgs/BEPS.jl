export Met, MetSeries, fill_met!

"""
# Fields
$(TYPEDFIELDS)
"""
@with_kw mutable struct Met <: AbstractFlux # 为套用相同结构
  "Inward shortwave radiation, `[W m⁻²]`"
  Rs::Cdouble = 0.0

  "(optional) Inward longwave radiation, `[W m⁻²]`"
  Rln_in::Cdouble = NaN

  "2m air temperature, `[°C]`"
  Tair::Cdouble = 0.0

  "Relative Humidity, `[%]`"
  RH::Cdouble = 0.0

  "precipitation, `[mm/h]`"
  Prcp::Cdouble = 0.0

  "Wind speed at measurement height z, `[m/s]`"
  Uz::Cdouble = 0.0
end
@DefFluxSeries MetSeries = Met

# Met(Rs, Rln_in, Tair, RH, Prcp, Uz) =
#   Met(; Rs, Rln_in, Tair, RH, Prcp, Uz)
function fill_met!(met::Met, forcing::MetSeries, i::Int)
    met.Rs = forcing.Rs[i]
    met.Tair = forcing.Tair[i]
    met.Prcp = forcing.Prcp[i]
    met.Uz = forcing.Uz[i]
    met.Rln_in = forcing.Rln_in[i]
    met.RH = forcing.RH[i]
end

# """
# - `Rs`: W m-2
# - `Tair`: degC
# - `Prcp`: mm
# - `Uz`: m/s
# - `RH`: relative humidity, %
# """
# function fill_met!(met::Met, Rs::FT, Tair::FT, Prcp::FT, Uz::FT, RH::FT)
#   met.Rs = Rs
#   met.Tair = Tair
#   met.Prcp = Prcp / 1000 # mm to m
#   met.Uz = Uz
#   met.Rln_in = NaN # use model longwave estimation unless overridden
#   met.RH = RH
# end

# function fill_met!(met::Met, d::DataFrame, k::Int=1; use_lrad::Bool=false)
#   RH = hasproperty(d, :RH) ? d.RH[k] : q2RH(d.qair[k], d.Tair[k])
#   fill_met!(met, d.Rs[k], d.Tair[k], d.Prcp[k], d.Uz[k], RH)
#   use_lrad && isfinite(d.Rln_in[k]) && (met.Rln_in = d.Rln_in[k])
# end

"""
    AirLayer{FT}

Struct to store environmental conditions in each air layer corresponds to one canopy layer.

# Fields
$(TYPEDFIELDS)

> Copied from Land.jl
"""
@with_kw mutable struct AirLayer{FT<:AbstractFloat}
  "Air temperature `[K]`"
  Tair::FT
  "Air density `[kg m⁻³]`"
  ρₐ::FT
  "Specific heat of air `[J kg⁻¹ K⁻¹]`"
  Cp_ca::FT
  "Vapor pressure deficit `[Pa]`"
  VPD::FT
  "Psychrometric constant `[Pa K⁻¹]`"
  γ::FT
  "Slope of Saturation vapor pressure es `[Pa K⁻¹]`"
  Δ::FT
  "Relative humility `[%]`"
  RH::FT
  "Wind speed `[m s⁻¹]`"
  wind::FT = FT(2)
end
