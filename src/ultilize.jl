using Dates, Statistics, DataFrames

parse_time(x::AbstractString) = DateTime(x, dateformat"yyyy-mm-ddTHH:MM:SSZ")


# 极高性能
@views function agg(A::AbstractArray, by; dims::Integer=3, fun=mean)
  @assert length(by) == size(A, dims)

  grps = sort!(unique(collect(by)))
  inds = Dict{eltype(grps),Vector{Int}}(g => Int[] for g in grps)
  @inbounds for i in eachindex(by)
    push!(inds[by[i]], i)
  end

  outsz = collect(size(A))
  outsz[dims] = length(grps)
  out = similar(A, Float64, Tuple(outsz))

  @inbounds for (j, g) in pairs(grps)
    r = fun(selectdim(A, dims, inds[g]); dims=dims)
    selectdim(out, dims, j) .= dropdims(r, dims=dims)
  end
  return out
end


## only for BEPS
function agg_daily(mat::AbstractMatrix, dates)
  dates_day = Date.(dates)
  r = agg(mat, dates_day; dims=1, fun=mean)
  r
end


function agg_daily(df_fluxes::AbstractDataFrame, dates)
  dates_day = Date.(dates)

  GPP_sim = df_fluxes.GPP
  ET_sim = df_fluxes.Trans + df_fluxes.Evap
  Hs_sim = df_fluxes.SH

  (;
    dates_day=unique(dates_day) |> sort,
    GPP_sim=agg(GPP_sim, dates_day; dims=1, fun=sum),
    ET_sim=agg(ET_sim, dates_day; dims=1, fun=sum),
    Hs_sim=agg(Hs_sim, dates_day; dims=1, fun=mean),
  )
end

export parse_time, agg, agg_daily
