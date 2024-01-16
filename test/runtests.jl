using Test
using BEPS
using BEPS: path_proj



include("test-clang.jl")
include("test-soil.jl")


@testset "besp_main julia" begin
  d = deserialize(path_proj("data/p1_meteo"))
  lai = readdlm(path_proj("examples/input/p1_lai.txt"))[:]

  par = (lon=120.5, lat=30.5, landcover=25, clumping=0.85,
    soil_type=8, Tsoil=2.2,
    soilwater=0.4115, snowdepth=0.0)

  @time df_jl, df_ET_jl = besp_main(d, lai, par; version="julia")
  # @time df_c, df_ET_c = besp_main(d, lai, par; version="c");
  r = sum(df_jl)

  @test abs(r.GPP - 2369.3039241523384) < 0.01
  @test abs(r.Evap - 748.8864673979658) < 0.01
  @test abs(r.Trans - 444.02624822679013) < 0.01

  # @test abs(r.GPP - 2369.3039943774525) < 1
  # @test abs(r.Evap - 748.8870728418663) < 1
  # @test abs(r.Trans - 444.0272797401601) < 1
end
