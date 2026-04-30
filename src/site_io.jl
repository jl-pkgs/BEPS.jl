export load_forcing, load_lai, load_site_info, load_multisite_data
export run_multisite   # extended overload: run_multisite(path; ...)


# ── Standard directory layout ────────────────────────────────────────────────
#
#  <data_dir>/
#  ├── site_info.csv          ← master metadata table  (required)
#  ├── <site_id>/
#  │   ├── forcing.txt        ← hourly meteorological forcing
#  │   └── lai.txt            ← daily LAI
#  └── ...
#
#  site_info.csv columns (required):
#    site_id  lon  lat  VegType  SoilType  clumping  Tsoil0  θ0  z_snow0
#  optional columns:
#    forcing_file  lai_file      ← override default paths per site
#    weight                      ← per-site weight for joint optimization
#
#  Forcing file (tab/space/comma separated, first row = header):
#    day  hour  Rs(or rad)  Tair(or tem)  RH(or hum/qair)  rain(or pre)  Uz(or wind)  [Rln_in]
#
#  LAI file (two accepted layouts):
#    a) one value per line: 365 or 366 rows
#    b) all values on a single row (space/tab separated)
# ─────────────────────────────────────────────────────────────────────────────


"""
    load_forcing(path) -> DataFrame

Load an hourly meteorological forcing file.

Accepted formats (auto-detected):
- Tab-/space-/comma-separated with a header row
- Required columns (after rename via `standardize_forcing_columns`):
  `day, hour, Rs/rad, Tair/tem, RH/hum/qair, rain/pre, Uz/wind`
- Optional: `Rln_in` / `lrad` (longwave radiation)

Returns a `DataFrame` with canonical column names ready for `beps_modern`.
"""
function load_forcing(path::AbstractString)::DataFrame
  isfile(path) || throw(ArgumentError("Forcing file not found: $path"))
  ext = lowercase(splitext(path)[2])
  delim = ext == ".csv" ? ',' : '\t'

  data, header = readdlm(path, delim, Any, header=true)
  cols = Symbol.(vec(header))
  # 尝试将每列转为数值
  d = DataFrame(
    [col => _coerce_numeric(data[:, i]) for (i, col) in enumerate(cols)]...
  )
  d = standardize_forcing_columns(d)
  # day/hour 必须为 Int（s_coszs 要求）
  d.day  = Int.(d.day)
  d.hour = Int.(d.hour)
  d
end

# 将 Any 数组尽量转为 Float64，失败时保留原始类型
function _coerce_numeric(col::AbstractVector)
  try
    return Float64.(col)
  catch
    return col
  end
end


"""
    load_lai(path) -> Vector{Float64}

Load a daily LAI file.

Accepted layouts:
- One value per line (365 or 366 rows)
- All values on a single space/tab-separated row

Returns a `Vector{Float64}`.
"""
function load_lai(path::AbstractString)::Vector{Float64}
  isfile(path) || throw(ArgumentError("LAI file not found: $path"))
  raw = readdlm(path)
  vec(Float64.(raw))
end


"""
    load_site_info(path) -> DataFrame

Load a site metadata CSV file.

Required columns: `site_id, lon, lat, VegType, SoilType, clumping, Tsoil0, θ0, z_snow0`

Optional columns:
- `forcing_file` : path to forcing file (absolute or relative to CSV directory)
- `lai_file`     : path to LAI file
- `weight`       : per-site weight for joint optimization (default 1.0)

Example `site_info.csv`:
```
site_id,lon,lat,VegType,SoilType,clumping,Tsoil0,θ0,z_snow0
CA-Obs,104.69,53.99,1,4,0.65,2.0,0.35,0.0
CA-Qfo,-74.34,49.69,1,3,0.60,1.5,0.40,0.0
```
"""
function load_site_info(path::AbstractString)::DataFrame
  isfile(path) || throw(ArgumentError("site_info file not found: $path"))
  data, header = readdlm(path, ',', header=true)
  cols = Symbol.(vec(header))
  df = DataFrame(data, cols)

  required = [:site_id, :lon, :lat, :VegType, :SoilType, :clumping, :Tsoil0, :θ0, :z_snow0]
  missing_cols = [string(c) for c in required if !(c in propertynames(df))]
  !isempty(missing_cols) &&
    throw(ArgumentError("site_info missing columns: $(join(missing_cols, ", "))"))

  # 强制类型转换
  df.site_id = string.(df.site_id)
  for col in [:lon, :lat, :clumping, :Tsoil0, :θ0, :z_snow0]
    df[!, col] = Float64.(df[!, col])
  end
  for col in [:VegType, :SoilType]
    df[!, col] = Int.(df[!, col])
  end

  # 可选列默认值
  :weight       in propertynames(df) || (df.weight = ones(nrow(df)))
  :forcing_file in propertynames(df) || (df.forcing_file = fill("", nrow(df)))
  :lai_file     in propertynames(df) || (df.lai_file     = fill("", nrow(df)))

  df
end


"""
    load_multisite_data(site_info; data_dir="") -> (forcing_dict, lai_dict)

Load forcing and LAI data for all sites described in `site_info`.

Data lookup order per site:
1. If `site_info.forcing_file` / `site_info.lai_file` is non-empty, use that path
   (relative paths resolved against `data_dir`).
2. Otherwise, look for `<data_dir>/<site_id>/forcing.txt` and
   `<data_dir>/<site_id>/lai.txt`.

Returns:
- `forcing_dict` : `Dict{String, DataFrame}`
- `lai_dict`     : `Dict{String, Vector{Float64}}`

Sites with missing files emit a warning and are skipped.
"""
function load_multisite_data(site_info::DataFrame;
    data_dir::AbstractString="")::Tuple{Dict{String,DataFrame},Dict{String,Vector{Float64}}}

  forcing_dict = Dict{String, DataFrame}()
  lai_dict     = Dict{String, Vector{Float64}}()

  for row in eachrow(site_info)
    id = string(row.site_id)

    # 解析强迫数据路径
    f_path = if !isempty(string(row.forcing_file))
      _resolve(string(row.forcing_file), data_dir)
    else
      joinpath(data_dir, id, "forcing.txt")
    end

    # 解析 LAI 路径
    l_path = if !isempty(string(row.lai_file))
      _resolve(string(row.lai_file), data_dir)
    else
      joinpath(data_dir, id, "lai.txt")
    end

    ok = true
    if !isfile(f_path)
      @warn "Forcing file not found for site $id: $f_path"
      ok = false
    end
    if !isfile(l_path)
      @warn "LAI file not found for site $id: $l_path"
      ok = false
    end
    ok || continue

    try
      forcing_dict[id] = load_forcing(f_path)
      lai_dict[id]     = load_lai(l_path)
    catch e
      @warn "Failed to load data for site $id" exception=e
    end
  end

  forcing_dict, lai_dict
end


"""
    run_multisite(site_info_path; data_dir, parallel, kw...) -> Dict

Convenience entry-point: read `site_info.csv` → load data → run all sites.

# Arguments
- `site_info_path` : path to `site_info.csv`
- `data_dir`       : base directory for resolving relative data paths
                     (defaults to the directory containing `site_info.csv`)
- `parallel`       : use `Threads.@threads` (default `true`)
- `kw`             : forwarded to `beps_modern`

# Returns
`Dict{String, NamedTuple}` where each value is `(; df_out, df_ET, Tsoil, θ)`.

# Typical workflow
```julia
using BEPS
results = run_multisite("path/to/site_info.csv")
```
"""
function run_multisite(site_info_path::AbstractString;
    data_dir::AbstractString=dirname(abspath(site_info_path)),
    parallel::Bool=true,
    kw...)

  site_info = load_site_info(site_info_path)
  forcing_dict, lai_dict = load_multisite_data(site_info; data_dir)
  run_multisite(site_info, forcing_dict, lai_dict; parallel, kw...)
end


# Resolve a path: if absolute, use as-is; else join with data_dir
_resolve(path, data_dir) = isabspath(path) ? path : joinpath(data_dir, path)
