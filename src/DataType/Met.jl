export Met

"""
# Fields
$(TYPEDFIELDS)
"""
@with_kw mutable struct Met
  "Inward shortwave radiation, `[W m⁻²]`"
  Rs::Cdouble = 0.0

  "Inward longwave radiation, `[W m⁻²]`"
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

# Met(Rs, Rln_in, Tair, RH, Prcp, Uz) =
#   Met(; Rs, Rln_in, Tair, RH, Prcp, Uz)

"""
# Arguments
- `Rs`: W m-2
- `Tair`: degC
- `Prcp`: mm
- `Uz`: m/s
- `RH`: relative humidity, %
"""
function fill_met!(met::Met,
  Rs::FT, Tair::FT, Prcp::FT, Uz::FT, RH::FT)

  met.Rs = Rs
  met.Tair = Tair
  met.Prcp = Prcp / 1000 # mm to m
  met.Uz = Uz
  met.Rln_in = NaN # use model longwave estimation unless overridden
  met.RH = RH
end

function fill_met!(met::Met, d::DataFrame, k::Int=1; use_lrad::Bool=false)
  RH = hasproperty(d, :RH) ? d.RH[k] : q2RH(d.qair[k], d.Tair[k])
  fill_met!(met, d.Rs[k], d.Tair[k], d.Prcp[k], d.Uz[k], RH)
  use_lrad && isfinite(d.Rln_in[k]) && (met.Rln_in = d.Rln_in[k])
end


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
  "Retive humility `[%]`"
  RH::FT
  "Wind speed `[m s⁻¹]`"
  wind::FT = FT(2)
end



const FORCING_RENAME_MAP = (
  :rad => :Rs,
  :tem => :Tair,
  :pre => :Prcp,
  :wind => :Uz,
  :hum => :qair,
  :lrad => :Rln_in,
)

function standardize_forcing_columns(d::DataFrame)
  cols = Set(Symbol.(names(d)))
  d2 = d

  for (old_col, new_col) in FORCING_RENAME_MAP
    if (old_col in cols) && !(new_col in cols)
      d2 === d && (d2 = copy(d))
      rename!(d2, old_col => new_col)
      push!(cols, new_col)
    end
  end

  required = (:day, :hour, :Rs, :Tair, :Prcp, :Uz)
  missing_cols = [string(c) for c in required if !(c in cols)]
  !isempty(missing_cols) &&
    throw(ArgumentError("Missing forcing columns: $(join(missing_cols, ", "))"))

  if !(:RH in cols) && !(:qair in cols)
    throw(ArgumentError("Missing humidity column: require RH or qair"))
  end

  if !(:Rln_in in cols)
    d2 === d && (d2 = copy(d))
    d2.Rln_in = fill(NaN, nrow(d2))
  end
  d2
end
