# BEPS
using RTableTools
using BEPS, DataFrames, Test, Dates

##
LAI = readdlm(path_proj("examples/input/p1_lai.txt"))[:]
forcing = deserialize(path_proj("data/p1_forcing"))
dates = DateTime(2010):Hour(1):DateTime(2010, 12, 31, 23)
ntime = length(dates)

VegType, SoilType = 25, 8
kw = (lon=120.5, lat=30.5, clumping=0.85, fix_snowpack=false)

ps = ParamBEPS(VegType, SoilType)
Ta = Float64(forcing.Tair[1])
state, _ = setup(ps; Ta, Tsoil=2.2, θ0=0.4115, z_snow=0.0)

df_jl, df_ET_jl, states_jl = beps_main(forcing, LAI, dates; ps, state, kw..., version="julia", verbose=false)


## 1. 测试C与Julia版本的结果是否一致
@time df_jl, df_ET_jl, states_jl = beps_main(forcing, LAI, dates; ps, state, kw..., version="julia")
@time df_c, df_ET_c, states_c = beps_main(forcing, LAI, dates; ps, state, kw..., version="c")
r = sum(df_jl)

df_diff = abs.(df_jl .- df_c)
df_diff_perc = abs.(df_jl .- df_c) ./ df_c .* 100

df_diff_perc = df_diff_perc[:, Cols(:GPP, :Evap, :Trans)]
l = maximum(df_diff_perc)

@show l
@show _nanmaximum(l)
@test true

# @test _nanmaximum(l) < 1.5 # SH, 1.48%的误差, current 0.09%
@test _nanmaximum(l) < 2.5 # GPP, 2.38%的误差, current 0.09%


## 2. 定位差异较大的时刻，输出数据进行对比分析
inds = findall(df_diff_perc.Trans .> 1)

if length(inds) > 0
  I = inds[1]
  index = max(1, I - 5):min(ntime, I + 5)
  @warn "debug: large difference"
  dat_c = cbind(forcing[index], df_c[index, 11:end])
  dat_jl = cbind(forcing[index], df_jl[index, 11:end])
  fwrite(dat_c, "dat_c.csv")
  fwrite(dat_jl, "dat_jl.csv")
end

## 3. 定位蒸发差异最大的时刻，输出数据进行对比分析
i_spike_evap = argmax(abs.(df_jl.Evap .- df_c.Evap))
@show i_spike_evap
@show df_jl.Evap[i_spike_evap]   # Julia绝对值
@show df_c.Evap[i_spike_evap]    # C绝对值
@show df_diff.Evap[i_spike_evap] # 绝对差


## 4. 绘图展示结果
using Plots
using Ipaper
gr(framestyle=:box)

function plot_var(var=:SH)
  n = size(df_jl, 1)
  inds = 1:n

  y_jl = df_jl[inds, var]
  y_c = df_c[inds, var]
  diff_perc = (y_jl .- y_c) ./ y_c .* 100
  plot(diff_perc; title="$var", label=nothing)
end

vars = names(df_jl)[[1:4; 8; 12:16]]
ps_plots = map(plot_var, vars)

plot(ps_plots..., layout=(3, 4), size=(1400, 800))
write_fig("Figure1_bias_of_julia-version-5deg.png", 15, 8; show=false)

## 土壤温度的变化
depths = [0.05, 0.10, 0.20, 0.40, 1.25] |> cumsum |> x -> round.(x, digits=4)
Tg = states_jl.vectors.Tsoil_c  # (5, ntime) matrix

p = plot(; size=(1400, 700), title="Soil temperature")
plot!(forcing.Tair; label="Tair")
for i = 1:5
  label = "depth: $(depths[i])"
  plot!(Tg[i, :]; label)
end
p
write_fig("Figure2_variation_of_Tg.png", 15, 8; show=false)
