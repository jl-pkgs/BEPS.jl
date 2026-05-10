## 1. 读取驱动数据
using BEPS, RTableTools, DataFrames, Dates, ModelParams
FT = Float64

indir = "z:/GitHub/jl-pkgs/ChinaFlux2026" |> path_mnt
st_full = fread("$indir/data/st_flux12.csv")

f = "$indir/data-raw/BEPS/Forcing_Met_Hourly_BEPS_Forest_sp12_hourly_v20260507.csv"
@time FORCING = fread("$indir/data-raw/BEPS/Forcing_Met_Hourly_BEPS_Forest_sp12_hourly_v20260507.csv")
replace_missing!(FORCING)
SITES = unique(FORCING.site)


##
SITE = SITES[8]
d_forcing = FORCING[FORCING.site.==SITE, 2:end]
rename!(d_forcing, [:Ta_canopy => :Tair, :RH_canopy => :RH, :WS_canopy => :Uz])
(; Tair, RH, Uz, Rs, Rln_in, Prcp) = d_forcing

dates_local = parse_time.(d_forcing.time)
dates_UTC = dates_local .- Hour(8) # [local] -> [UTC]
ntime = length(dates_local)

forcing = MetSeries(; ntime, Rs, Rln_in, Tair, RH, Uz, Prcp)

# 率定所需数据
f = "$indir/data-raw/BEPS/Fluxes/$(SITE)_FluxLAISoil_daily_v20260510.csv" |> path_mnt
FluxALL = fread(f)
replace_missing!(FluxALL)
rename!(FluxALL, :LAI_glass_G005 => :lai, :GPP => :GPP_obs, :ET => :ET_obs, :Hs => :Hs_obs)
(; lai) = FluxALL


## 2. 初始化模型参数和状态变量
st = st_full[findfirst(st_full.site .== SITE), :]
(; lon, lat, VegType, SoilType, z_Uz, z_overstory) = st
depths_SM = map(x -> parse(Int, x), split(st.z_SM, ",")) ./ 100 # 
depths_TS = map(x -> parse(Int, x), split(st.z_TS, ",")) ./ 100 # 

model = ParamBEPS(VegType, SoilType)
model.veg.z_wind = z_Uz
model.veg.z_canopy_o = z_overstory

model.veg = ParamVeg{FT}(;
  VCmax25=56.25,
  g1_w=4.8,
  Ω=0.58,
)

state = InitState0(model, forcing)
@time df_fluxes, df_ET, states, caches = simulate(forcing, lai, dates_UTC;
  ps=model, state, lon, lat);
@time gof, data_sim, data_obs = BEPS_GOF(df_fluxes, states, dates_local, FluxALL;
  depths_SM, depths_TS)
gof

# fwrite(cbind(; time=dates, df_ET = df_ET .* 3600), "./Project_ChinaFlux/df_ET.csv")
## 加入参数优化模块
opts = [
  (; path=[:r_drainage], name=:r_drainage, value=0.5),
  (; path=[:veg, :Ω], name=:Ω, value=0.58),
  (; path=[:veg, :g1_w], name=:g1_w, value=4.8),
  (; path=[:veg, :g0_w], name=:g0_w, value=0.0175),
  (; path=[:veg, :VCmax25], name=:VCmax25, value=56.25),
] |> DataFrame
paths = opts.path
kw_loss = (; lon, lat, depths_SM, depths_TS, FluxDay=FluxALL)

@time theta_opt = optim(model, forcing, lai, dates_UTC; paths, maxn=1000, kw_loss...)
opts.theta_opt = theta_opt
opts

##
include("main_vis.jl")
plot_result(data_sim, data_obs)
