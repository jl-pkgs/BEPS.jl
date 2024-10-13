@testset "rainfall_stage1_jl" begin
  Tair = 5.0  # Air temperature in Celsius
  prcp = 10.0/1000  # Precipitation in mm
  perc_water = Layer2{Float64}(0.1, 0.2)  # Percentage water
  m_water = Layer2{Float64}(0.01, 0.02)  # Current water mass
  m_water_pre = Layer2{Float64}(0.05, 0.005)  # Previous water mass
  lai_o = 1.0  # Leaf area index of overstory
  lai_u = 1.0  # Leaf area index of understory
  Ω = 0.5  # Clumping index

  # Call the function with the test inputs
  prcp_g = rainfall_stage1_jl(Tair, prcp, perc_water, m_water, m_water_pre, lai_o, lai_u, Ω)

  r2 = clang.rainfall_stage1(Tair, prcp, perc_water, m_water, m_water_pre, lai_o, lai_u, Ω)

  @test prcp_g == r2
end
