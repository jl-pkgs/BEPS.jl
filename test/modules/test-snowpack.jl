using BEPS, Test
using BEPS: snowpack_stage3

# @testset "snowpack_stage1" 
begin
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
  m_snow = Layer3(0.1, 0.2, 0.3)    # by reference
  perc_snow = Layer3(0.1, 0.2, 0.3) # by reference
  area_snow = Layer2(0.1, 0.2)      # by reference
  @test call_jl() == call_c()

  Tair = 1.0
  m_snow = Layer3(0.1, 0.2, 0.3)
  perc_snow = Layer3(0.1, 0.2, 0.3)
  area_snow = Layer2(0.1, 0.2)
  r1 = call_jl()
  r2 = call_c()

  @test call_jl() == call_c()

  Tair = 0.0001 # near zero, 精度不足存在一定问题
  m_snow = Layer3(0.1, 0.2, 0.3)   # value was changed
  perc_snow = Layer3(0.1, 0.2, 0.3)
  area_snow = Layer2(0.1, 0.2)
  @test call_jl() == call_c()
end

@testset "snowpack_stage3" begin
  Ω = 0.5
  depth_snow = 0.5
  Tair = -1.0

  ρ_snow = 700.0        # [kg m-3]
  depth_snow = 0.06     # [m]
  depth_water = 0.02    # [m]

  # 融化
  Tsnow_last = -1.0
  Tsnow = 3.0

  for depth_snow = [0.01, 0.03, 0.06]
    r_c = snowpack_stage3(Tair, Tsnow, Tsnow_last, ρ_snow, depth_snow, depth_water, Layer3(0.1, 0.2, 3.0))
    r_jl = snowpack_stage3_jl(Tair, Tsnow, Tsnow_last, ρ_snow, depth_snow, depth_water, Layer3(0.1, 0.2, 3.0))
    @test r_c == r_jl
  end
  
  ## 结冻
  Tsnow_last = 3.0
  Tsnow = -1.0
  r_c = snowpack_stage3(Tair, Tsnow, Tsnow_last, ρ_snow, depth_snow, depth_water, Layer3(0.1, 0.2, 3.0))
  r_jl = snowpack_stage3_jl(Tair, Tsnow, Tsnow_last, ρ_snow, depth_snow, depth_water, Layer3(0.1, 0.2, 3.0))
  #! C结冻存在错误
  @test r_jl == (0.06000080057281437, 0.019999453267346284)
end
