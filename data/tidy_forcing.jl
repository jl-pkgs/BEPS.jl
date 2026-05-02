using RTableTools, Dates

dates = DateTime(2010):Hour(1):DateTime(2010, 12, 31, 23)
dates_daily = Date(2010):Day(1):Date(2010, 12, 31)

d_forcing = fread("./examples/input/p1_meteo.txt")

(; Rs, Tair, q, Prcp, Uz) = d_forcing
Tair = Tair .- 5.0

ntime = length(dates)
RH = q2RH.(q, Tair)
Rln_in = fill(NaN, ntime)
forcing = MetSeries{Float64}(; ntime, Rs, Tair, RH, Prcp=Prcp / 1000, Uz, Rln_in)
serialize(path_proj("data/p1_forcing"), forcing)
