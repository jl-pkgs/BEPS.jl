export run_multisite, beps_multisite_optimize


"""
    run_multisite(sites, forcing_dict, lai_dict; parallel=true, kw...) -> Dict

Run BEPS for multiple sites in parallel (via `Threads.@threads`) or serial.

# Arguments
- `sites`        : `DataFrame` with columns `site_id, lon, lat, VegType, SoilType,
                   clumping, Tsoil0, θ0, z_snow0`
- `forcing_dict` : `Dict{String, DataFrame}` — hourly meteorological forcing per site
- `lai_dict`     : `Dict{String, Vector}` — daily LAI per site
- `parallel`     : use multi-threading (default `true`)
- `kw`           : extra keyword arguments forwarded to `beps_modern`

# Returns
`Dict{String, NamedTuple}` where each value is `(; df_out, df_ET, Tsoil, θ)`.
"""
function run_multisite(
    sites::DataFrame,
    forcing_dict::Dict,
    lai_dict::Dict;
    parallel::Bool=true,
    kw...)

  n  = nrow(sites)
  out = Dict{String,NamedTuple}()
  lk  = ReentrantLock()

  function _run_one(i)
    row = sites[i, :]
    id  = string(row.site_id)
    haskey(forcing_dict, id) || (@warn "No forcing for site $id"; return)
    haskey(lai_dict, id)     || (@warn "No LAI for site $id";     return)

    df_out, df_ET, Tsoil, θ = beps_modern(
      forcing_dict[id], lai_dict[id];
      lon      = Float64(row.lon),
      lat      = Float64(row.lat),
      VegType  = Int(row.VegType),
      SoilType = Int(row.SoilType),
      clumping = Float64(row.clumping),
      Tsoil0   = Float64(row.Tsoil0),
      θ0       = Float64(row.θ0),
      z_snow0  = Float64(row.z_snow0),
      kw...   # let caller control verbose and other options
    )
    lock(lk) do
      out[id] = (; df_out, df_ET, Tsoil, θ)
    end
  end

  if parallel
    Threads.@threads for i in 1:n
      _run_one(i)
    end
  else
    for i in 1:n
      _run_one(i)
    end
  end
  out
end


"""
    beps_multisite_optimize(sites, forcing_dict, lai_dict, obs_dict, model;
                            col_sim, weights, normalize, opt_paths, kw...)

Multi-site joint optimization of **shared** parameters in `model`.

All sites use the same `model` (shared parameter set, e.g. same vegetation type).
Per-site initial conditions come from the `sites` DataFrame.
The objective is the weighted average of per-site normalized RMSE.

# Arguments
- `sites`        : `DataFrame` (columns: `site_id, lon, lat, VegType, SoilType,
                   clumping, Tsoil0, θ0, z_snow0`)
- `forcing_dict` : `Dict{String, DataFrame}` — hourly forcing per site
- `lai_dict`     : `Dict{String, Vector}` — daily LAI per site
- `obs_dict`     : observations per site; either
                   `Dict{String, AbstractVector}` (single column), or
                   `Dict{String, Dict{Symbol, AbstractVector}}` (multi-objective per site)
- `model`        : `ParamBEPS` — shared parameter container, modified in-place
- `col_sim`      : output column when `obs_dict` values are plain vectors (default `:GPP`)
- `weights`      : `Dict{String, Float64}` — per-site weights (default equal)
- `normalize`    : normalize RMSE by `std(obs)` for each site (default `true`)
- `opt_paths`    : restrict which parameter paths are optimized
"""
function beps_multisite_optimize(
    sites::DataFrame,
    forcing_dict::Dict,
    lai_dict::Dict,
    obs_dict::Dict,
    model::ParamBEPS;
    col_sim::Symbol=:GPP,
    weights::Dict{String,Float64}=Dict(string(r.site_id) => 1.0 / nrow(sites)
                                       for r in eachrow(sites)),
    normalize::Bool=true,
    opt_paths=nothing,
    maxn=2000, kstop=5, f_reltol=0.0001, x_reltol=0.0001, seed=1,
    n_complex=5, size_complex=nothing, size_simplex=nothing,
    n_evolu=nothing, n_pop=nothing, verbose=false,
    kwargs...)

  x0, lb, ub, paths = _filter_opt_info(model, opt_paths)

  function cal_func(x)
    m = deepcopy(model)
    update!(m, paths, x)

    J = 0.0;  W = 0.0
    for row in eachrow(sites)
      id = string(row.site_id)
      (haskey(forcing_dict, id) && haskey(lai_dict, id) && haskey(obs_dict, id)) || continue
      w = get(weights, id, 1.0)
      w == 0 && continue

      df_out = df_ET = nothing
      try
        df_out, df_ET, _, _ = beps_modern(
          forcing_dict[id], lai_dict[id];
          model    = m,
          lon      = Float64(row.lon),
          lat      = Float64(row.lat),
          VegType  = Int(row.VegType),
          SoilType = Int(row.SoilType),
          clumping = Float64(row.clumping),
          Tsoil0   = Float64(row.Tsoil0),
          θ0       = Float64(row.θ0),
          z_snow0  = Float64(row.z_snow0),
          kwargs...  # let caller control verbose and other options
        )
      catch e
        e isa DomainError || @warn "beps_multisite_optimize: site $id error" exception=e
        J += w * 1e6;  W += w
        continue
      end

      obs_site = obs_dict[id]
      site_obj = if obs_site isa Dict
        # multi-objective per site
        isempty(obs_site) && throw(ArgumentError("obs_dict[\"$id\"] must not be empty"))
        sim_dict = Dict(col => _get_col(df_out, df_ET, col) for col in keys(obs_site))
        wts_site = Dict(col => 1.0 / length(obs_site) for col in keys(obs_site))
        _weighted_rmse(obs_site, sim_dict, wts_site, normalize)
      else
        sim  = _get_col(df_out, df_ET, col_sim)
        n    = min(length(obs_site), length(sim))
        rmse = of_RMSE(obs_site[1:n], sim[1:n])
        if normalize
          valid_obs = filter(!isnan, obs_site[1:n])
          length(valid_obs) < 2 && (J += w * 1e6; W += w; continue)
          σ = std(valid_obs)
          (!isfinite(σ) || σ == 0) ? Inf : rmse / σ
        else
          rmse
        end
      end

      J += w * site_obj;  W += w
    end
    W == 0 ? Inf : J / W
  end

  # 多站点目标函数较重，禁用函数内并行以防嵌套线程
  bestx, bestf, _ = sceua(cal_func, x0, lb, ub;
    _sceua_kw(; maxn, kstop, f_reltol, x_reltol, seed, n_complex, verbose,
               parallel=false, size_complex, size_simplex, n_evolu, n_pop)...)
  update!(model, paths, bestx)
  model, bestf
end
