using BEPS, RTableTools, DataFrames, Dates, ModelParams
using Ipaper


function agg_daily(mat::AbstractMatrix, dates)
  dates_day = Date.(dates)
  r = apply(mat, 2; by=dates_day, fun=nanmean)
  collect(r')
end

function agg_daily(df_fluxes::AbstractDataFrame, dates)
  dates_day = Date.(dates)

  GPP_sim = df_fluxes.GPP
  ET_sim = df_fluxes.Trans + df_fluxes.Evap
  Hs_sim = df_fluxes.SH

  fun = nansum
  (;
    GPP_sim=apply(GPP_sim, 1; by=dates_day, fun),
    ET_sim=apply(ET_sim, 1; by=dates_day, fun),
    Hs_sim=apply(Hs_sim, 1; by=dates_day, fun=nanmean),
    dates_day=unique_sort(dates_day)
  )
end

parse_time(x::AbstractString) = DateTime(x, dateformat"yyyy-mm-ddTHH:MM:SSZ")
