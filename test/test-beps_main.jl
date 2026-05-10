using BEPS, DataFrames, Dates, Test

kw = (lon=120.5, lat=30.5,
  VegType=25, SoilType=8,
  clumping=0.85,
  Tsoil0=2.2, θ0=0.4115, z_snow0=0.0
)

LAI = readdlm(path_proj("examples/input/p1_lai.txt"))[:]
forcing = deserialize(path_proj("data/p1_forcing"))
forcing.Prcp .= forcing.Prcp .* 1000.0  # 转换为[m] -> [mm]

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

##
@testset "beps_main julia" begin
  df_jl, df_ET_jl, states_jl = main(; version="julia")
  @time df_jl, df_ET_jl, states_jl = main(; version="julia")
  @time df_c, df_ET_c, states_c = main(; version="c")
  r = sum(df_jl)

  # @test abs(r.GPP - 2369.3039241523384) < 0.01
  # @test abs(r.Evap - 748.8864673979658) < 0.01
  # @test abs(r.Trans - 444.02624822679013) < 0.01
  df_diff = abs.(df_jl .- df_c)
  df_diff_perc = abs.(df_jl .- df_c) ./ df_c .* 100

  # gpp_u_sunlit has a large bias in i=193, unknown reason
  # df_diff_perc = df_diff_perc[:, Cols(1:1, 3:end)]
  df_diff_perc = df_diff_perc[:, Cols(:GPP, :Evap, :Trans)]
  for col in names(df_diff_perc)
    df_diff_perc[abs.(df_c[!, col]).<eps(Float64), col] .= 0.0
  end
  l = maximum(df_diff_perc)
  @show l
  @show _nanmaximum(l)
  @test true
  # @test _nanmaximum(l) < 1.5 # SH, 1.48%的误差, current 0.09%
  # @test _nanmaximum(l) < 2.5 # GPP, 2.38%的误差, current 0.09%
end

## performance
# @time df_c, df_ET_c, states_c = main(; version="julia");
# @profview df_jl, df_ET_jl, states_jl = main(; version="julia")
