export of_NSE, of_KGE, of_R2, of_Bias, of_MAE
export evaluate_site, evaluate_multisite

using Statistics: mean, std, cor


"""
    of_NSE(obs, sim) -> Float64

Nash-Sutcliffe Efficiency. Range (-∞, 1]. Perfect fit = 1.
"""
function of_NSE(obs::AbstractVector, sim::AbstractVector)
  valid = @. !isnan(obs) & !isnan(sim)
  o, s = obs[valid], sim[valid]
  isempty(o) && return NaN
  denom = sum((o .- mean(o)) .^ 2)
  denom == 0 && return NaN
  1.0 - sum((o .- s) .^ 2) / denom
end


"""
    of_KGE(obs, sim) -> Float64

Kling-Gupta Efficiency (Gupta et al. 2009). Range (-∞, 1]. Perfect = 1.
KGE = 1 - sqrt((r-1)² + (α-1)² + (β-1)²), where α = σ_s/σ_o, β = μ_s/μ_o.
"""
function of_KGE(obs::AbstractVector, sim::AbstractVector)
  valid = @. !isnan(obs) & !isnan(sim)
  o, s = obs[valid], sim[valid]
  length(o) < 2 && return NaN
  μo, μs = mean(o), mean(s)
  σo, σs = std(o), std(s)
  (σo == 0 || μo == 0) && return NaN
  r = cor(o, s)
  !isfinite(r) && return NaN
  α = σs / σo  # variability ratio
  β = μs / μo  # bias ratio
  1.0 - sqrt((r - 1)^2 + (α - 1)^2 + (β - 1)^2)
end


"""
    of_R2(obs, sim) -> Float64

Coefficient of determination (R²). Range [0, 1].
"""
function of_R2(obs::AbstractVector, sim::AbstractVector)
  valid = @. !isnan(obs) & !isnan(sim)
  o, s = obs[valid], sim[valid]
  length(o) < 2 && return NaN
  (std(o) == 0 || std(s) == 0) && return NaN
  r = cor(o, s)
  !isfinite(r) && return NaN
  r^2
end


"""
    of_Bias(obs, sim) -> Float64

Mean relative bias (%). Positive = overestimate, negative = underestimate.
Uses `mean(obs)` as denominator; returns `NaN` when `mean(obs) == 0`.
"""
function of_Bias(obs::AbstractVector, sim::AbstractVector)
  valid = @. !isnan(obs) & !isnan(sim)
  o, s = obs[valid], sim[valid]
  isempty(o) && return NaN
  μo = mean(o)
  μo == 0 && return NaN
  (mean(s) - μo) / abs(μo) * 100.0
end


"""
    of_MAE(obs, sim) -> Float64

Mean absolute error.
"""
function of_MAE(obs::AbstractVector, sim::AbstractVector)
  valid = @. !isnan(obs) & !isnan(sim)
  o, s = obs[valid], sim[valid]
  isempty(o) && return NaN
  mean(abs.(o .- s))
end


"""
    evaluate_site(obs, sim) -> NamedTuple

Compute performance metrics for one site/variable pair.
Returns `(; RMSE, NSE, KGE, R2, Bias, MAE)`.
"""
function evaluate_site(obs::AbstractVector, sim::AbstractVector)
  (;
    RMSE = of_RMSE(obs, sim),
    NSE  = of_NSE(obs, sim),
    KGE  = of_KGE(obs, sim),
    R2   = of_R2(obs, sim),
    Bias = of_Bias(obs, sim),
    MAE  = of_MAE(obs, sim),
  )
end


"""
    evaluate_multisite(results, obs_dict; col=:GPP) -> DataFrame

Evaluate model performance across multiple sites.

# Arguments
- `results`  : `Dict{String, NamedTuple}` returned by `run_multisite`
- `obs_dict` : `Dict{String, AbstractVector}` — observed values per site
- `col`      : output column to compare (default `:GPP`); searched first in `df_out`,
               then in `df_ET`

# Returns
A `DataFrame` with columns `site, RMSE, NSE, KGE, R2, Bias, MAE` (one row per site).
"""
function evaluate_multisite(
    results::Dict,
    obs_dict::Dict;
    col::Symbol=:GPP,
)
  rows = NamedTuple[]
  for site_id in sort(collect(keys(results)))
    haskey(obs_dict, site_id) || continue
    res = results[site_id]
    obs = obs_dict[site_id]
    df_out = res.df_out
    df_ET  = res.df_ET

    sim = if col in propertynames(df_out)
      df_out[!, col]
    elseif col in propertynames(df_ET)
      df_ET[!, col]
    else
      @warn "Column $col not found for site $site_id, skipping"
      continue
    end

    n = min(length(obs), length(sim))
    length(obs) != length(sim) &&
      @warn "obs/sim length mismatch for site $site_id ($(length(obs)) vs $(length(sim))), using first $n"
    metrics = evaluate_site(obs[1:n], sim[1:n])
    push!(rows, (; site=site_id, metrics...))
  end
  DataFrame(rows)
end
