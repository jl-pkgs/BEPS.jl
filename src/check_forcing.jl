export sanitize_forcing!, normalize_flux_obs!


function _interp_invalid!(x::AbstractVector{<:AbstractFloat})
  good = findall(isfinite, x)
  isempty(good) && return x

  first_good = first(good)
  last_good = last(good)
  x[begin:first_good-1] .= x[first_good]
  x[last_good+1:end] .= x[last_good]

  good = findall(isfinite, x)
  @inbounds for k in 1:length(good)-1
    i, j = good[k], good[k+1]
    j == i + 1 && continue
    xi, xj = x[i], x[j]
    for t in i+1:j-1
      x[t] = xi + (xj - xi) * (t - i) / (j - i)
    end
  end
  return x
end

function _sanitize_column!(d::AbstractDataFrame, col::Symbol;
  lo=-Inf, hi=Inf, fill=NaN, interp=true, clamp_lo=-Inf, clamp_hi=Inf)

  col in propertynames(d) || return 0

  x = d[!, col]
  bad = (.!isfinite.(x)) .| (x .< lo) .| (x .> hi)
  x[bad] .= fill
  interp && _interp_invalid!(x)

  if isfinite(clamp_lo) && isfinite(clamp_hi)
    x .= clamp.(x, clamp_lo, clamp_hi)
  elseif isfinite(clamp_lo)
    x .= max.(x, clamp_lo)
  elseif isfinite(clamp_hi)
    x .= min.(x, clamp_hi)
  end
  return count(bad)
end

function sanitize_forcing!(d::AbstractDataFrame)
  bad_Tair = _sanitize_column!(d, :Tair; lo=-50.0, hi=50.0)
  bad_RH = _sanitize_column!(d, :RH; lo=0.0, hi=100.0, clamp_lo=0.0, clamp_hi=100.0)
  bad_Uz = _sanitize_column!(d, :Uz; lo=0.0, clamp_lo=0.01)
  bad_Rs = _sanitize_column!(d, :Rs; lo=0.0, clamp_lo=0.0)
  bad_Prcp = _sanitize_column!(d, :Prcp; lo=0.0, fill=0.0, interp=false, clamp_lo=0.0)
  return (; bad_Tair, bad_RH, bad_Uz, bad_Rs, bad_Prcp)
end

function normalize_flux_obs!(d::AbstractDataFrame)
  names_d = propertynames(d)
  if :GPP_obs in names_d && mean(skipmissing(d.GPP_obs)) < 0
    d.GPP_obs .= -d.GPP_obs
  end
  return d
end

function _daily_obs_dates(d::AbstractDataFrame)
  :date in propertynames(d) || return nothing
  Date.(d.date)
end

function _align_daily_data(data_sim::DataFrame, data_obs::DataFrame)
  common_dates = intersect(data_sim.date, data_obs.date)
  sort!(common_dates)
  isempty(common_dates) && error("No overlapping dates between simulated and observed daily data")

  sim_idx = indexin(common_dates, data_sim.date)
  obs_idx = indexin(common_dates, data_obs.date)
  return data_sim[sim_idx, :], data_obs[obs_idx, :]
end
