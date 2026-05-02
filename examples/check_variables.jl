# BEPS
using RTableTools, Ipaper
using BEPS, DataFrames, Test, Dates, Plots
gr(framestyle=:box)

function build_forcing(; Tair_offset=0.0)
  d_forcing = fread("./examples/input/p1_meteo.txt")
  (; Rs, Tair, q, Prcp, Uz) = d_forcing
  Tair = Tair .+ Tair_offset

  ntime = length(dates)
  RH = q2RH.(q, Tair)
  Rln_in = fill(NaN, ntime)
  MetSeries{Float64}(; ntime, Rs, Tair, RH, Prcp=Prcp / 1000, Uz, Rln_in)
end


## 1. 跑模型
LAI = readdlm(path_proj("examples/input/p1_LAI.txt"))[:]
dates = DateTime(2010):Hour(1):DateTime(2010, 12, 31, 23)
ntime = length(dates)

VegType, SoilType = 25, 8
kw = (lon=120.5, lat=30.5, clumping=0.85, 
  fix_snowpack=true, fix_Ta_annual=true, 
  version="julia", verbose=false)

# forcing = deserialize(path_proj("data/p1_forcing"))
# forcing = build_forcing(; Tair_offset=-5.0)
forcing = build_forcing(; Tair_offset=0.0)

ps = ParamBEPS(VegType, SoilType)
Ta = Float64(forcing.Tair[1])
state, _ = setup(ps; Ta, Tsoil=2.2, θ0=0.4115, z_snow=0.0)

df_jl, df_ET_jl, states_jl = beps_main(forcing, LAI, dates; ps, state, kw...)

# 土壤温度的变化
depths = [0.05, 0.10, 0.20, 0.40, 1.25] |> cumsum |> x -> round.(x, digits=4)
Tsoil = states_jl.vectors.Tsoil_c  # (5, ntime) matrix


## Figure_Tsoil
p = plot(; size=(1400, 700), title="Soil temperature")
plot!(forcing.Tair; label="Tair")

for i = 1:5
  label = "depth: $(depths[i])"
  plot!(Tsoil[i, :]; label)
end
p
write_fig("Figure2_variation_of_Tsoil.png", 15, 8; show=false)


## check_snow
plot(
  plot(dates, states_jl.scalars.z_water .* 1000; title="Water depth (mm)", label="z_water"),
  plot(dates, states_jl.scalars.z_snow .* 1000; title="Snow depth (mm)", label="z_snow"),
  plot(dates, states_jl.scalars.ρ_snow; title="ρ_snow (kg m^3)", label="ρ_snow"),
  plot(dates, forcing.Prcp .* 1000; title="Precipitation (mm)", label="Prcp"),
  plot(dates, forcing.Tair; title="Air temperature (°C)", label="Tair"),
  size = (1400, 700),
  layout = (2, 3),
)
write_fig("Figure3_snow_process.png", 15, 8; show=false)


## check_SM
p = plot(; size=(1400, 700), title="Soil Moisture")
# plot!(forcing.Tair; label="Tair")
SM = states_jl.vectors.θ  # (5, ntime) matrix
for i = 1:5
  label = "depth: $(depths[i])"
  plot!(SM[i, :]; label)
end
p
write_fig("Figure2_variation_of_SM.png", 15, 8; show=false)
