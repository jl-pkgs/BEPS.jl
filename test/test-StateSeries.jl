using BEPS, Test

const SF = (:z_water, :z_snow, :r_rain_g, :f_soilwater) |> Val
const VF = (:θ, :Tsoil_c, :ETi, :G) |> Val
const NLAYER = 5

@testset "StateSeries" begin
    ntime = 8760
    out = StateSeries(SF, VF, NLAYER, ntime)
    display(out)

    VegType::Int = 25
    SoilType::Int = 8
    model = ParamBEPS(VegType, SoilType)
    state = setup(model)[1]

    t = 1
    save_state!(out, state, t, SF, VF)
    s = out[1]
    @test s[:θ] == [0.24, 0.3, 0.315, 0.33, 0.345]
end

## 最后把数据保存为nc
