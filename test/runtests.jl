using Test
using BEPS
using BEPS: path_proj


function nanmax(x)
  x = collect(x)
  x = x[.!isnan.(x)]
  maximum(x)
end

@testset "besp_main julia" begin
  d = deserialize(path_proj("data/p1_meteo"))
  lai = readdlm(path_proj("examples/input/p1_lai.txt"))[:]

  par = (lon=120.5, lat=30.5, landcover=25, clumping=0.85,
    soil_type=8, Tsoil=2.2,
    soilwater=0.4115, snowdepth=0.0)

  @time df_jl, df_ET_jl = besp_main(d, lai, par; version="julia")
  @time df_c, df_ET_c = besp_main(d, lai, par; version="c");
  r = sum(df_jl)

  @test abs(r.GPP - 2369.3039241523384) < 0.01
  @test abs(r.Evap - 748.8864673979658) < 0.01
  @test abs(r.Trans - 444.02624822679013) < 0.01

  df_diff = abs.(df_jl .- df_c)
  df_diff_perc = abs.(df_jl .- df_c) ./ df_c .* 100
  # max(df_diff_perc)
  
  l = max(df_diff_perc)
  @test nanmax(l) < 1.5 # SH, 1.48%的误差, current 0.09%
end


include("test-clang.jl")
include("test-soil.jl")
