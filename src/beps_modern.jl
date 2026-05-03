"""
    beps_modern(forcing, lai, dates; ps, state, ...)

Run BEPS simulation using the modern `ParamBEPS + StateBEPS` API.

# Arguments
- `forcing`     : `MetSeries` with hourly meteorological data
- `lai`         : daily LAI vector (length = number of days)
- `dates`       : `DateTime` vector matching `forcing` length
- `ps`          : `ParamBEPS` model parameters
- `state`       : `StateBEPS` initial state (copied internally, external value unchanged)
- `lon`, `lat`  : longitude and latitude [°]
- `clumping`    : canopy clumping index (default 0.85)
- `sm_obs`      : observed soil moisture [m³/m³], `nlayer × ntime` matrix; when provided,
                  `UpdateSoilMoisture` is skipped and each hourly θ is prescribed from data
- `VARS_EXPORT` : state variables to save (default `DEFAULT_VARS_EXPORT`)

# Returns
`(df_flux, df_ET, states)` — hourly flux DataFrames and a `StateSeries`
"""
function beps_modern(forcing::MetSeries, lai::Vector, dates::AbstractVector;
  ps::ParamBEPS, state::StateBEPS,
  lon::FT=120.0, lat::FT=20.0, clumping::FT=0.85,
  fix_snowpack=true, fix_annual_Ta=true,
  sm_obs::Union{Nothing, AbstractMatrix}=nothing,
  VARS_EXPORT::Vector{Symbol}=DEFAULT_VARS_EXPORT, kw...) where {FT<:AbstractFloat}

  state = deepcopy(state) # 避免外部状态被修改
  met = Met()
  mid_flux = Flux()
  mid_ET = ETFlux()
  cache = LeafCache()

  ntime = forcing.ntime
  fluxes = FluxSeries(; ntime)
  fluxes_ET = ETSeries(; ntime)

  SF, VF = split_vars(VARS_EXPORT)
  states = StateSeries(SF, VF, layer, ntime)

  Ta_annual = mean(forcing.Tair)
  jdays = dayofyear.(dates)
  hours = hour.(dates)

  fix_sm = sm_obs !== nothing

  for i = 1:ntime
    jday = jdays[i]
    hour = hours[i]

    if fix_sm
      state.θ_prev .= state.θ
      state.θ[1:5] .= @view sm_obs[:, i]
    end

    fill_met!(met, forcing, i) # 驱动数据
    k = ceil(Int, i / 24)
    _lai = lai[k]

    inter_prg_jl(jday, hour, lon, lat, _lai, clumping,
      met, ps, state, mid_flux, mid_ET, cache; fix_snowpack, fix_annual_Ta, Ta_annual, fix_sm)

    fluxes[i] = mid_flux
    fluxes_ET[i] = mid_ET
    save_state!(states, state, i, SF, VF)
  end
  DataFrame(fluxes), DataFrame(fluxes_ET), states
end
