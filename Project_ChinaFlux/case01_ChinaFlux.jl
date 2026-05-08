# 快为第一准则，先把模型跑起来
using BEPS
using RTableTools, DataFrames, Dates, ModelParams
using Ipaper


function agg_daily(df_fluxes, dates)
  dates_day = Date.(dates)

  GPP_sim = df_fluxes.GPP
  ET_sim = df_fluxes.Trans + df_fluxes.Evap
  fun = nansum
  (;
    GPP_sim=apply(GPP_sim, 1; by=dates_day, fun),
    ET_sim=apply(ET_sim, 1; by=dates_day, fun),
    dates_day=unique_sort(dates_day)
  )
end

parse_time(x::AbstractString) = DateTime(x, dateformat"yyyy-mm-ddTHH:MM:SSZ")

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
ntime = length(dates)

forcing = MetSeries(; ntime, Rs, Rln_in, Tair, RH, Uz, Prcp)

d_flux = FluxLAI[FluxLAI.site.==SITE, [:GPP, :ET, :LAI_glass_G005]]
rename!(d_flux, :LAI_glass_G005 => :lai, :GPP => :GPP_obs, :ET => :ET_obs)
(; lai, GPP_obs, ET_obs) = d_flux

GPP_obs = -GPP_obs # gC m-2 day-1, [GEE] -> [GPP]
ET_obs = ET_obs * 1e+06 / 86400 # [MJ m-2] -> [W m-2], mm


## 2. 初始化模型参数和状态变量
# 初始化模型参数
# VegType=1: ENF (Evergreen Needleleaf Forest)，对应千烟洲人工针叶林
VegType::Int = 1
SoilType::Int = 8
model = ParamBEPS(VegType, SoilType)

Ta = 20.0
Tsoil0 = Ta
θ0 = model.hydraulic.θ_vfc[1]  # 初始化为田间持水量（避免低于凋萎含水量）
z_snow0 = 0.0
state, model = setup(model; Ta, Tsoil=Float64(Tsoil0), θ0=Float64(θ0), z_snow=Float64(z_snow0))
state

@time df_fluxes, df_ET, states, caches = beps_modern(forcing, lai, dates; ps=model, state,
  lon=115.06, lat=26.74, clumping=0.62);  # ENF clumping=0.62
# 0.452331 seconds (1.14 M allocations: 53.704 MiB, 11.91% gc time)
# > 模型顺利开跑，0.45s跑完2年hourly 17544步长模拟
(GPP_sim, ET_sim, dates_day) = agg_daily(df_fluxes, dates)

fwrite(cbind(; time=dates, df_ET = df_ET .* 3600), "./Project_ChinaFlux/df_ET.csv")

##
gof_gpp = GOF(GPP_obs, GPP_sim)
gof_et = GOF(ET_obs, ET_sim)

DataFrame([
  (; var="GPP", gof_gpp...),
  (; var="ET", gof_et...)
])


##
using Plots

p1 = plot(dates_day, GPP_obs, label="GPP_obs")
plot!(p1, dates_day, GPP_sim, label="GPP_sim")

p2 = plot(dates_day, ET_obs, label="ET_obs")
plot!(p2, dates_day, ET_sim, label="ET_sim")

p = plot(p1, p2)
savefig(p, "a.png")
