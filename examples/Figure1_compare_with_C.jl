using BEPS, Plots, Dates
gr(; framestyle=:box)

##
forcing = deserialize(path_proj("data/p1_forcing"))
LAI = readdlm(path_proj("examples/input/p1_LAI.txt"))[:]
dates = DateTime(2010):Hour(1):DateTime(2010, 12, 31, 23)

VegType, SoilType = 25, 8
kw = (lon=120.5, lat=30.5, clumping=0.85, fix_snowpack=false)

ps = ParamBEPS(VegType, SoilType)
Ta = Float64(forcing.Tair[1])
state0, _ = setup(ps; Ta, Tsoil=2.2, θ0=0.4115, z_snow=0.0)


## Compare with C
@time df_jl, df_ET_jl, states_jl = besp_main(forcing, LAI, dates; ps, state=state0, kw..., version="julia")
@time df_c, df_ET_c, states_c = besp_main(forcing, LAI, dates; ps, state=state0, kw..., version="c")

sum(df_jl)
sum(df_c)

df_diff = abs.(df_jl .- df_c)
df_diff_perc = abs.(df_jl .- df_c) ./ df_c .* 100
maximum(df_diff)
maximum(df_diff_perc) # SH, 1.48%的误差


## Plot
function plot_var(var=:SH)
  n = size(df_jl, 1)
  inds = 1:n

  y_jl = df_jl[inds, var]
  y_c = df_c[inds, var]
  diff_perc = (y_jl .- y_c) ./ y_c .* 100

  plot(diff_perc; title="$var", label="bias (%)")
end

vars = names(df_jl)[[1:4; 8; 12:16]]
ps_plots = map(plot_var, vars)
plot(ps_plots..., layout=(3, 4), size=(1400, 800))
# savefig("images/Figure1_bias_of_julia-version.png")

# # using ProfileView
# # @profview_allocs df_jl, df_ET_jl, _ = besp_main(forcing, LAI, dates; ps, state=state0, kw..., version="julia");
# # @profview df_jl, df_ET_jl, _ = besp_main(forcing, LAI, dates; ps, state=state0, kw..., version="julia");
