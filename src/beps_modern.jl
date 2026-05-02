"""
    beps_modern(d, lai; model, VegType, SoilType, ...)

Run BEPS simulation using the modern `ParamBEPS + StateBEPS` API.

# Arguments
- `d`   : hourly forcing `DataFrame`. Required columns:

  | Column   | Unit    | Description                        |
  |----------|---------|------------------------------------|
  | `day`    | —       | Day of year (1–366)                |
  | `hour`   | —       | Hour of day (1–24)                 |
  | `Rs`     | W m⁻²   | Downward shortwave radiation       |
  | `Tair`   | °C      | Air temperature                    |
  | `Prcp`   | mm h⁻¹  | Precipitation                      |
  | `Uz`     | m s⁻¹   | Wind speed                         |
  | `qair`   | g kg⁻¹  | Specific humidity *(or `RH` [%])*  |
  | `Rln_in` | W m⁻²   | Downward longwave radiation *(optional)* |

- `lai`     : daily LAI vector (length = number of days)
- `model`   : `ParamBEPS` object; if `nothing`, built from `VegType`/`SoilType`
- `VegType` : IGBP vegetation type (used when `model=nothing`)
- `SoilType`: USDA soil texture class (used when `model=nothing`)

# Returns
`(df_out, df_ET, output_Tsoil, output_θ)`
"""
function beps_modern(d::DataFrame, lai::Vector;
  model::Union{Nothing,ParamBEPS}=nothing,
  lon::FT=120.0, lat::FT=20.0,
  VegType::Int=25, SoilType::Int=8, clumping::FT=0.85,
  Tsoil0::FT=2.2, θ0::FT=0.4115, z_snow0::FT=0.0,
  r_drainage::FT=0.5, r_root_decay::FT=0.95,
  fix_snowpack=true, verbose=true, kw...) where {FT<:AbstractFloat}

  Ta = Float64(d.Tair[1])

  # 构建 ParamBEPS 和 StateBEPS
  if isnothing(model)
    # 与 besp_main 使用相同的初始化路径，确保结果一致
    _, state, ps = setup_model(VegType, SoilType;
      version="julia", Ta, Tsoil=Float64(Tsoil0), θ0=Float64(θ0),
      z_snow=Float64(z_snow0), r_drainage, r_root_decay, FT)
  else
    ps = model
    state, _ = setup(ps; Ta, Tsoil=Float64(Tsoil0), θ0=Float64(θ0), z_snow=Float64(z_snow0))
  end

  met = Met()
  mid_res = Flux()
  mid_ET = ETFlux()
  Ra = Radiation()
  cache = LeafCache()

  ntime = size(d, 1)
  vars = fieldnames(Flux) |> collect
  vars_ET = fieldnames(ETFlux) |> collect

  df_out = DataFrame(zeros(ntime, length(vars)), vars)
  df_ET = DataFrame(zeros(ntime, length(vars_ET)), vars_ET)
  output_Tsoil = zeros(ntime, layer)
  output_θ = zeros(ntime, layer)

  for i = 1:ntime
    jday = d.day[i]
    hour = d.hour[i]
    CosZs = s_coszs(jday, hour, lat, lon)

    _day = ceil(Int, i / 24)
    (mod(_day, 50) == 0 && hour == 1 && verbose) && println("Day = $_day")

    _lai = lai[_day]
    fill_met!(met, d, i)

    inter_prg_jl(jday, hour, lon, lat, _lai, clumping, Ra, met, ps, state,
      mid_res, mid_ET, cache; fix_snowpack)

    output_Tsoil[i, :] .= state.Tsoil_c[1:layer]
    output_θ[i, :] .= state.θ[1:layer]
    fill_res!(df_out, mid_res, i)
    fill_res!(df_ET, mid_ET, i)
  end

  df_out, df_ET, output_Tsoil, output_θ
end

# ETSeries
# FluxSeries
# StateSeries
