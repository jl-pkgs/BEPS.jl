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

function sanitize_forcing!(d::AbstractDataFrame)
  stats = (; bad_Tair=0, bad_RH=0, bad_Uz=0, bad_Rs=0, bad_Prcp=0)

  if :Tair in propertynames(d)
    bad = (.!isfinite.(d.Tair)) .| (d.Tair .< -50.0) .| (d.Tair .> 50.0)
    d.Tair[bad] .= NaN
    _interp_invalid!(d.Tair)
    stats = merge(stats, (; bad_Tair=count(bad)))
  end
  if :RH in propertynames(d)
    bad = (.!isfinite.(d.RH)) .| (d.RH .< 0.0) .| (d.RH .> 100.0)
    d.RH[bad] .= NaN
    _interp_invalid!(d.RH)
    d.RH .= clamp.(d.RH, 0.0, 100.0)
    stats = merge(stats, (; bad_RH=count(bad)))
  end
  if :Uz in propertynames(d)
    bad = (.!isfinite.(d.Uz)) .| (d.Uz .< 0.0)
    d.Uz[bad] .= NaN
    _interp_invalid!(d.Uz)
    d.Uz .= max.(d.Uz, 0.01)
    stats = merge(stats, (; bad_Uz=count(bad)))
  end
  if :Rs in propertynames(d)
    bad = (.!isfinite.(d.Rs)) .| (d.Rs .< 0.0)
    d.Rs[bad] .= NaN
    _interp_invalid!(d.Rs)
    d.Rs .= max.(d.Rs, 0.0)
    stats = merge(stats, (; bad_Rs=count(bad)))
  end
  if :Prcp in propertynames(d)
    bad = (.!isfinite.(d.Prcp)) .| (d.Prcp .< 0.0)
    d.Prcp[bad] .= 0.0
    d.Prcp .= max.(d.Prcp, 0.0)
    stats = merge(stats, (; bad_Prcp=count(bad)))
  end
  return stats
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
