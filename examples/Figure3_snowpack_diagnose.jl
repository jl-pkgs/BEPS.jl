using BEPS, Plots, Dates

d = deserialize("data/p1_meteo")
d.tem = d.tem .- 5.0 # 为了测试融雪模块
lai = readdlm("examples/input/p1_lai.txt")[:]

par = (lon=120.5, lat=30.5, landcover=25, clumping=0.85,
  soil_type=8, Tsoil=2.2,
  soilwater=0.4115, snowdepth=0.0)

function plot_var(var=:SH; data=df_out)
  dates = DateTime(2010):Hour(1):DateTime(2010, 12, 31, 23)
  ticks = DateTime(2010):Month(2):DateTime(2010, 12, 31, 23)
  ticklabels = Dates.format.(ticks, "mm")
  xticks = ticks, ticklabels

  n = size(data, 1)
  inds = 1:n
  # inds = 1:24*90
  y_jl = data[inds, var]
  plot(dates, y_jl; title="$var", xticks)
end

function plot_Figure3(fout)
  gr(framestyle=:box)
  p_tem = plot_var("tem"; data=d)
  hline!(p_tem, [0], color="red")
  plot(
    plot_var("z_water"),
    plot_var("z_snow"),
    plot_var("ρ_snow"), 
    plot_var("pre"; data=d),
    p_tem,
    size=(800, 600)
  )
  savefig(fout)
end

df_out, df_ET_jl, Tsoil, θ = besp_main(d, lai, par;version="julia")
plot_Figure3("docs/images/Figure3_BEPS_snowpack_Julia_v0.1.7.png")

df_out, df_ET_jl, Tsoil, θ = besp_main(d, lai, par; version="julia", fix_snowpack=false)
plot_Figure3("docs/images/Figure3_BEPS_snowpack_BEPS_v4.10.png")
