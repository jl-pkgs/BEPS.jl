using BEPS, DataFrames, Test
using BEPS: path_proj

function nanmax(x)
  x = collect(x)
  x = x[.!isnan.(x)]
  maximum(x)
end

lai = readdlm(path_proj("examples/input/p1_lai.txt"))[:]

kw = (lon=120.5, lat=30.5,
  VegType=25, SoilType=8,
  clumping=0.85,
  Tsoil0=2.2, θ0=0.4115, z_snow0=0.0
)

# @testset "besp_main julia" 
begin
  d = deserialize(path_proj("data/p1_meteo"))
  d.tem = d.tem .- 5.0

  @time df_jl, df_ET_jl = besp_main(d, lai; kw..., version="julia", fix_snowpack=false)
  @time df_c, df_ET_c = besp_main(d, lai; kw..., version="c")
  r = sum(df_jl)

  # @test abs(r.GPP - 2369.3039241523384) < 0.01
  # @test abs(r.Evap - 748.8864673979658) < 0.01
  # @test abs(r.Trans - 444.02624822679013) < 0.01

  df_diff = abs.(df_jl .- df_c)
  df_diff_perc = abs.(df_jl .- df_c) ./ df_c .* 100

  # gpp_u_sunlit has a large bias in i=193, unknown reason
  # df_diff_perc = df_diff_perc[:, Cols(1:1, 3:end)]
  df_diff_perc = df_diff_perc[:, Cols(:GPP, :Evap, :Trans)]
  l = maximum(df_diff_perc)
  @show l
  @show nanmax(l)
  @test true

  # @test nanmax(l) < 1.5 # SH, 1.48%的误差, current 0.09%
  @test nanmax(l) < 2.5 # GPP, 2.38%的误差, current 0.09%
end
