# 快为第一准则，先把模型跑起来
include("main_pkgs.jl")

## 1. 读取驱动数据
indir = "Z:/GitHub/jl-pkgs/ChinaFlux2026/data-raw/BEPS"
indir = "/mnt/z/GitHub/jl-pkgs/ChinaFlux2026/data-raw/BEPS"

@time FORCING = fread("$indir/Forcing_Met_Hourly_BEPS_Forest_sp12_hourly_v20260507.csv")
replace_missing!(FORCING)
SITES = unique(FORCING.site)

@time FluxLAI = fread("$indir/Forcing_FluxLAI_Daily_BEPS_Forest_sp12_v20260506.csv")
replace_missing!(FluxLAI)

##
SITE = SITES[8]
d_forcing = FORCING[FORCING.site.==SITE, 2:end]
rename!(d_forcing, [:Ta_canopy => :Tair, :RH_canopy => :RH, :WS_canopy => :Uz])
(; Tair, RH, Uz, Rs, Rln_in, Prcp) = d_forcing
# Prcp = Prcp ./ 1000.0  # 转换为[mm] -> [m]

dates = parse_time.(d_forcing.time)
dates_model = dates .- Hour(8) # [local] -> [UTC]
ntime = length(dates)

forcing = MetSeries(; ntime, Rs, Rln_in, Tair, RH, Uz, Prcp)

d_flux = FluxLAI[FluxLAI.site.==SITE, [:GPP, :ET, :LAI_glass_G005, :Hs]]
rename!(d_flux, :LAI_glass_G005 => :lai, :GPP => :GPP_obs, :ET => :ET_obs, :Hs => :Hs_obs)
(; lai, GPP_obs, ET_obs, Hs_obs) = d_flux

GPP_obs = -GPP_obs # gC m-2 day-1, [GEE] -> [GPP]
ET_obs = ET_obs * 1e+06 / 86400 # [MJ m-2] -> [W m-2], ET的最终单位是mm
Hs_obs = Hs_obs * 1e+06 / 86400 # [MJ m-2] -> [W m-2], Hs的最终单位是W m-2

## 2. 初始化模型参数和状态变量
# 初始化模型参数
# VegType=1: ENF (Evergreen Needleleaf Forest)，对应千烟洲人工针叶林
VegType::Int = 1
SoilType::Int = 8
model = ParamBEPS(VegType, SoilType)
model.veg.z_wind = 39.6
model.veg.VCmax25 = 56.25
model.veg.g1_w = 4.8
clumping = 0.58

Ta = forcing.Tair[1]
Tsoil0 = Ta
θ0 = model.hydraulic.θ_vfc[1]  # 初始化为田间持水量（避免低于凋萎含水量）
z_snow0 = 0.0
state, model = setup(model; Ta, Tsoil=Float64(Tsoil0), θ0=Float64(θ0), z_snow=Float64(z_snow0))
state

@time df_fluxes, df_ET, states, caches = simulate(forcing, lai, dates_model; ps=model, state,
  lon=115.06, lat=26.74, clumping);
# 0.452331 seconds (1.14 M allocations: 53.704 MiB, 11.91% gc time)
# > 模型顺利开跑，0.45s跑完2年hourly 17544步长模拟
(GPP_sim, ET_sim, Hs_sim, dates_day) = agg_daily(df_fluxes, dates)
# fwrite(cbind(; time=dates, df_ET = df_ET .* 3600), "./Project_ChinaFlux/df_ET.csv")

gof = DataFrame([
  (; var="GPP", GOF(GPP_obs, GPP_sim)...),
  (; var="ET", GOF(ET_obs, ET_sim)...),
  (; var="Hs", GOF(Hs_obs, Hs_sim)...)
])


##
using Plots
gr(; framestyle=:box)

p1 = plot(dates_day, GPP_obs, label="GPP_obs")
plot!(p1, dates_day, GPP_sim, label="GPP_sim")

p2 = plot(dates_day, ET_obs, label="ET_obs")
plot!(p2, dates_day, ET_sim, label="ET_sim")

p3 = plot(dates_day, Hs_obs, label="Hs_obs")
plot!(p3, dates_day, Hs_sim, label="Hs_sim")

p = plot(p1, p2, p3, layout=(3, 1), size=(1400, 900))
savefig(p, "Figure1_fluxes.png")

## SM
labels = map(i -> "L$i", 1:5)
SM_day = agg_daily(states.vectors.θ, dates)

p = plot(dates_day, SM_day; label=labels, size = (1400, 400))
savefig(p, "Figure1_SM.png")

TS_day = agg_daily(states.vectors.Tsoil_c, dates)
p = plot(dates_day, TS_day; label=labels, size=(1400, 400))
savefig(p, "Figure1_TS.png")
