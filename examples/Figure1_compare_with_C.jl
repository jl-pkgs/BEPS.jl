using BEPS, Plots, Dates

d = deserialize("data/p1_meteo")
d.tem = d.tem .- 5.0 # 为了测试融雪模块
lai = readdlm("examples/input/p1_lai.txt")[:]

par = (lon=120.5, lat=30.5, landcover=25, clumping=0.85,
  soil_type=8,
  Tsoil=2.2,
  soilwater=0.4115,
  snowdepth=0.0)


## Compare with C
@time df_jl, df_ET_jl = besp_main(d, lai, par; version="julia", fix_snowpack=false);
@time df_c, df_ET_c = besp_main(d, lai, par; version="c");

sum(df_jl)
sum(df_c)

df_diff = abs.(df_jl .- df_c)
df_diff_perc = abs.(df_jl .- df_c) ./ df_c .* 100
max(df_diff)
max(df_diff_perc) # SH, 1.48%的误差


## Plot
function plot_var(var=:SH)
  n = size(df_jl, 1)
  inds = 1:n

  y_jl = df_jl[inds, var]
  y_c = df_c[inds, var]
  diff_perc = (y_jl .- y_c) ./ y_c .* 100

  plot(diff_perc; title="$var", label="bias (%)") # 
end

vars = names(df_jl)[[1:4; 8; 12:16]]
ps = map(plot_var, vars)
plot(ps..., layout=(3, 4), size=(1400, 800))
# savefig("images/Figure1_bias_of_julia-version.png")

# # using ProfileView
# # @time df_jl, df_ET_jl = besp_main(d, lai, par; version="julia");
# # @profview df_jl, df_ET_jl = besp_main(d, lai, par; version="julia");
# @profview_allocs df_jl, df_ET_jl = besp_main(d, lai, par; version="julia");
# @profview df_jl, df_ET_jl = besp_main(d, lai, par; version="julia");
