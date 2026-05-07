# 快为第一准则，先把模型跑起来
using BEPS
using RTableTools, DataFrames, Dates

parse_time(x::AbstractString) = DateTime(x, dateformat"yyyy-mm-ddTHH:MM:SSZ")

## 1. 读取驱动数据

@time FORCING = fread("Z:/GitHub/jl-pkgs/ChinaFlux2026/data-raw/BEPS/Forcing_Met_Hourly_BEPS_Forest_sp12_hourly_v20260507.csv")
@time LAI = fread("Z:/GitHub/jl-pkgs/ChinaFlux2026/data-raw/BEPS/Forcing_FluxLAI_Daily_BEPS_Forest_sp12_v20260506.csv")

replace_missing!(FORCING)
replace_missing!(LAI)
SITES = unique(FORCING.site)

##
SITE = SITES[1]
d_forcing = FORCING[FORCING.site .== SITE, 2:end]
rename!(d_forcing, [:Ta_canopy => :Tair, :RH_canopy => :RH, :WS_canopy => :Uz])
(; Tair, RH, Uz, Rs, Rln_in, Prcp) = d_forcing
dates = parse_time.(d_forcing.time)
ntime = length(dates)

forcing = MetSeries(; ntime, Rs, Rln_in, Tair, RH, Uz, Prcp)
lai = LAI[LAI.site.==SITE, :LAI_glass_G005]


## 2. 初始化模型参数和状态变量
# 初始化模型参数
VegType::Int = 25
SoilType::Int = 8

model = ParamBEPS(VegType, SoilType)

## 初始化状态变量
Ta = 20.0
Tsoil0 = Ta
θ0 = model.hydraulic.θ_vfc[1] .* 0.6
z_snow0 = 0.0
state, model = setup(model; Ta, Tsoil=Float64(Tsoil0), θ0=Float64(θ0), z_snow=Float64(z_snow0))
state

## 开干
@time df_fluxes, df_ET, states, caches = beps_modern(forcing, lai, dates; ps=model, state);
# 0.452331 seconds (1.14 M allocations: 53.704 MiB, 11.91% gc time)
# > 模型顺利开跑，0.45s跑完2年hourly 17544步长模拟
