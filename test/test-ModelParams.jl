using BEPS, Test


@testset "Model Parameters" begin
  # xs = ParamSoilHydraulicLayers{Float64,4}()
  model = ParamBEPS{Float64}(; N=5)
  params = parameters(model)
  display(model)

  paths = [
    [:r_drainage],
    [:hydraulic, :b, 4]
  ]
  values = [0.4, 4.0]

  update!(model, paths, values; params)
  @test model.r_drainage == 0.4
  @test model.:hydraulic.b[4] == 4.0
end

# @testset "ParamSoilHydraulicLayers" begin
#   x = ParamSoilHydraulic{Float64}()
#   xs = ParamSoilHydraulicLayers{Float64,4}()
#   length(get_bounds(x)) == 6
#   length(get_bounds(xs)) == 24
# end
