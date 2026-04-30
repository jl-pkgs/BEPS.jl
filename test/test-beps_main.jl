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

@testset "besp_main julia" begin
  d = deserialize(path_proj("data/p1_meteo"))
  d.tem = d.tem .- 5.0

  df_jl, df_ET_jl = besp_main(d, lai; kw..., version="julia", verbose=false, fix_snowpack=false)
  @time df_jl, df_ET_jl = besp_main(d, lai; kw..., version="julia", fix_snowpack=false);
  @time df_c, df_ET_c = besp_main(d, lai; kw..., version="c")
  r = sum(df_jl)

  # @test abs(r.GPP - 2369.3039241523384) < 0.01
  # @test abs(r.Evap - 748.8864673979658) < 0.01
  # @test abs(r.Trans - 444.02624822679013) < 0.01

  df_diff_perc = abs.(df_jl .- df_c) ./ (abs.(df_c) .+ 1e-6) .* 100
  df_diff_perc = df_diff_perc[:, Cols(:GPP, :Evap, :Trans)]
  l = maximum(df_diff_perc)
  @show nanmax(l)
  @test nanmax(l) < 2.5  # Julia V1 should match C within 2.5%
end
