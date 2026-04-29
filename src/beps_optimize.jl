"""
    beps_optimize(d, lai, model, obs; col_sim, opt_paths, kw...) -> (model, best_rmse)

Single-site, single-objective optimization (RMSE).

# Arguments
- `obs`       : observed vector
- `col_sim`   : output column to compare (default `:ET`)
- `opt_paths` : restrict optimization to a subset of parameter paths, e.g.
                `[[:veg, :VCmax25], [:veg, :g1_w]]`; `nothing` = all params.
"""
function beps_optimize(d::DataFrame, lai::Vector, model::ParamBEPS,
    obs::AbstractVector;
    col_sim::Symbol=:ET,
    opt_paths=nothing,
    # SCE-UA
    maxn=2000, kstop=5, f_reltol=0.0001, x_reltol=0.0001, seed=1,
    n_complex=5, size_complex=nothing, size_simplex=nothing,
    n_evolu=nothing, n_pop=nothing, verbose=false, parallel=true,
    # beps_modern
    kwargs...)

  x0, lb, ub, paths = _filter_opt_info(model, opt_paths)

  function cal_func(x)
    m = deepcopy(model)
    update!(m, paths, x)
    df_out = df_ET = nothing
    try
      df_out, df_ET, _, _ = beps_modern(d, lai; model=m, kwargs...)
    catch e
      e isa DomainError || @warn "beps_optimize: unexpected error" exception=e
      return Inf
    end
    sim = _get_col(df_out, df_ET, col_sim)
    of_RMSE(obs, sim)
  end

  bestx, bestf, _ = sceua(cal_func, x0, lb, ub;
    _sceua_kw(; maxn, kstop, f_reltol, x_reltol, seed, n_complex, verbose, parallel,
               size_complex, size_simplex, n_evolu, n_pop)...)
  update!(model, paths, bestx)
  model, bestf
end


"""
    beps_optimize(d, lai, model, obs::Dict; weights, normalize, opt_paths, kw...)

Multi-objective single-site optimization.

# Arguments
- `obs`       : `Dict{Symbol, AbstractVector}` — e.g. `Dict(:GPP => gpp_obs, :LH => lh_obs)`
- `weights`   : per-variable weights (default equal); must sum to > 0
- `normalize` : divide each RMSE by `std(obs[col])` to remove unit dependence (default `true`)
- `opt_paths` : restrict optimization to a subset of parameter paths
"""
function beps_optimize(d::DataFrame, lai::Vector, model::ParamBEPS,
    obs::Dict{Symbol,<:AbstractVector};
    weights::Dict{Symbol,Float64}=Dict(k => 1.0 / length(obs) for k in keys(obs)),
    normalize::Bool=true,
    opt_paths=nothing,
    maxn=2000, kstop=5, f_reltol=0.0001, x_reltol=0.0001, seed=1,
    n_complex=5, size_complex=nothing, size_simplex=nothing,
    n_evolu=nothing, n_pop=nothing, verbose=false, parallel=true,
    kwargs...)

  x0, lb, ub, paths = _filter_opt_info(model, opt_paths)

  function cal_func(x)
    m = deepcopy(model)
    update!(m, paths, x)
    df_out = df_ET = nothing
    try
      df_out, df_ET, _, _ = beps_modern(d, lai; model=m, kwargs...)
    catch e
      e isa DomainError || @warn "beps_optimize (multi-obj): unexpected error" exception=e
      return Inf
    end
    sim_dict = Dict(col => _get_col(df_out, df_ET, col) for col in keys(obs))
    _weighted_rmse(obs, sim_dict, weights, normalize)
  end

  bestx, bestf, _ = sceua(cal_func, x0, lb, ub;
    _sceua_kw(; maxn, kstop, f_reltol, x_reltol, seed, n_complex, verbose, parallel,
               size_complex, size_simplex, n_evolu, n_pop)...)
  update!(model, paths, bestx)
  model, bestf
end


# ── Shared helpers (used by both optimize and multisite) ────────────────────────

# Filter get_opt_info result to a subset of paths
function _filter_opt_info(model, opt_paths)
  x0_all, lb_all, ub_all, paths_all = get_opt_info(model)
  isnothing(opt_paths) && return x0_all, lb_all, ub_all, paths_all
  inds = findall(p -> p in opt_paths, paths_all)
  isempty(inds) && error("opt_paths matched no parameters in model")
  x0_all[inds], lb_all[inds], ub_all[inds], paths_all[inds]
end

# Extract a column from df_out or df_ET (error if not found)
function _get_col(df_out, df_ET, col::Symbol)
  col in propertynames(df_out) && return df_out[!, col]
  col in propertynames(df_ET)  && return df_ET[!, col]
  error("Column $col not found in output")
end

# Weighted normalized RMSE across multiple variables
function _weighted_rmse(obs_dict, sim_dict, weights, normalize)
  J = 0.0; W = 0.0
  for (col, obs) in obs_dict
    w   = get(weights, col, 1.0)
    sim = sim_dict[col]
    n   = min(length(obs), length(sim))
    rmse = of_RMSE(obs[1:n], sim[1:n])
    σ   = normalize ? std(filter(!isnan, obs[1:n])) : 1.0
    σ == 0 && continue
    J += w * rmse / σ;  W += w
  end
  W == 0 ? Inf : J / W
end

# Build SCE-UA keyword NamedTuple (drops nothing-valued extras)
function _sceua_kw(; maxn, kstop, f_reltol, x_reltol, seed, n_complex,
                     verbose, parallel, size_complex, size_simplex, n_evolu, n_pop)
  base = (; maxn, kstop, f_reltol, x_reltol, seed, n_complex, verbose, parallel)
  ext  = (; size_complex, size_simplex, n_evolu, n_pop)
  merge(base, NamedTuple(k => v for (k, v) in pairs(ext) if !isnothing(v)))
end
