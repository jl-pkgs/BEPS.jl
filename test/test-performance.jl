using BEPS, DataFrames, Dates, Test

kw = (lon=120.5, lat=30.5,
  VegType=25, SoilType=8,
  clumping=0.85,
  Tsoil0=2.2, θ0=0.4115, z_snow0=0.0
)

LAI = readdlm(path_proj("examples/input/p1_lai.txt"))[:]
forcing = deserialize(path_proj("data/p1_forcing"))
dates = DateTime(2010):Hour(1):DateTime(2010, 12, 31, 23)

function main(; version="julia")
  beps_main(forcing, LAI, dates; kw..., version, verbose=false, fix_snowpack=false, fix_Ta_annual=false)
end

## tidy forcing
@testset "beps_main " begin
  @time df_jl, df_ET_jl, states_jl = main(; version="julia")
  r = sum(df_jl)
  @test isapprox(r.GPP, 2146.110; atol=2)
  @test isapprox(r.Evap, 62.5378; atol=1)
end

## performance
# @time df_jl, df_ET_jl, states_jl = main(; version="julia");
@profview for i in 1:10
  df_jl, df_ET_jl, states_jl = main(; version="julia")
end
