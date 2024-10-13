using BEPS, Test

@testset "snowpack_stage1" begin
  function call_jl()
    ρ_snow = Ref{Float64}(0.3)
    albedo_v_snow = Ref{Float64}(0.8)
    albedo_n_snow = Ref{Float64}(0.7)

    r = snowpack_stage1_jl(Tair, prcp, lai_o, lai_u, Ω, m_snow_pre, m_snow, perc_snow, area_snow, depth_snow, ρ_snow, albedo_v_snow, albedo_n_snow)

    r, ρ_snow[], albedo_v_snow[], albedo_n_snow[]
  end
  
  function call_c()
    ρ_snow = Ref{Float64}(0.3)
    albedo_v_snow = Ref{Float64}(0.8)
    albedo_n_snow = Ref{Float64}(0.7)

    r = clang.snowpack_stage1(Tair, prcp, lai_o, lai_u, Ω, m_snow_pre, m_snow, perc_snow, area_snow, depth_snow, ρ_snow, albedo_v_snow, albedo_n_snow)

    r, ρ_snow[], albedo_v_snow[], albedo_n_snow[]
  end

  prcp = 10.0
  lai_o = 1.0
  lai_u = 1.0
  Ω = 0.5
  m_snow_pre = Layer3(0.1, 0.2, 0.3)
  depth_snow = 0.5

  Tair = -1.0
  m_snow = Layer3(0.1, 0.2, 0.3)   # value was changed
  perc_snow = Layer3(0.1, 0.2, 0.3)
  area_snow = Layer2(0.1, 0.2)
  @test call_jl() == call_c()

  Tair = 1.0
  m_snow = Layer3(0.1, 0.2, 0.3)   # value was changed
  perc_snow = Layer3(0.1, 0.2, 0.3)
  area_snow = Layer2(0.1, 0.2)
  @test call_jl() == call_c()

  Tair = 0.0001 # near zero, 精度不足存在一定问题
  m_snow = Layer3(0.1, 0.2, 0.3)   # value was changed
  perc_snow = Layer3(0.1, 0.2, 0.3)
  area_snow = Layer2(0.1, 0.2)
  r1 = call_jl()
  r2 = call_c()
  @test call_jl() == call_c()
end

# julia>   m_snow
# Layer3{Float64}
# o: Float64 1.4519019656603828e6
# u: Float64 880623.1964169953
# g: Float64 1.357475437922622e6
