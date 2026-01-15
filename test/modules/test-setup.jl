using Test
using BEPS
import BEPS: setup, setup_jl, setup_c,
  Init_Soil_Parameters, Init_Soil_T_θ!, UpdateSoilMoisture, UpdateRootFraction!


@testset "setup" begin

  @testset "setup(VegType, SoilType; ...)" begin
    st, ps = setup(25, 8; Ta=20.0, θ0=0.3)

    @test st.n_layer == 5
    @test ps.N == 5
    @test st.dz[1:5] == ps.dz

    # 检查温湿度初始化
    @test all(st.θ[1:5] .> 0)
    @test all(st.Tsoil_c[1:5] .≈ 20.0)

    # 检查根系分布
    @test sum(st.f_root[1:5]) ≈ 1.0 atol = 0.01
  end


  @testset "setup(ps::ParamBEPS; ...)" begin
    ps = ParamBEPS(6, 3)  # DBF, sandy loam
    st, ps_ret = setup(ps; Ta=15.0, θ0=0.25)

    @test ps_ret === ps  # 返回相同的参数对象
    @test st.n_layer == ps.N
    @test all(st.Tsoil_c[1:5] .≈ 15.0)
    @test st.θ[2] ≈ 0.25  # 第2层应该等于 θ0
  end


  @testset "setup(soil::Soil)" begin
    soil = Soil()
    Init_Soil_Parameters(soil, 25, 8, 0.95)
    Init_Soil_T_θ!(soil, 18.0, 20.0, 0.35, 0.0)

    st, ps = setup(soil)

    # 状态变量应该与原 soil 一致
    @test st.n_layer == soil.n_layer
    @test st.θ[1:5] == soil.θ[1:5]
    @test st.Tsoil_c[1:5] == soil.Tsoil_c[1:5]
    @test st.ice_ratio[1:5] == soil.ice_ratio[1:5]

    # 参数应该与原 soil 一致
    @test ps.N == Int(soil.n_layer)
    @test ps.hydraulic.θ_sat == soil.θ_sat[1:5]
    @test ps.hydraulic.K_sat == soil.K_sat[1:5]
  end


  @testset "setup with keyword arguments" begin
    params = (VegType=25, SoilType=8, Ta=20.0, θ0=0.3)
    st, ps = setup(; params...)

    @test st.n_layer == 5
    @test ps.N == 5
  end


  @testset "UpdateSoilMoisture compatibility" begin
    st, ps = setup(25, 8; Ta=20.0, θ0=0.3)
    θ_before = copy(st.θ[1:5])

    # 添加降水
    st.r_rain_g = 0.0001  # [m/s]
    UpdateSoilMoisture(st, ps, 360.0)

    # 土壤湿度应该有变化
    @test st.θ[1:5] != θ_before
    @test all(st.θ[1:5] .> 0)
  end


  @testset "DBF/EBF parameters" begin
    # 注意：is_bforest 需要在 InitParam_Veg 中正确设置
    # 当前 ParamBEPS 构造函数使用 veg.is_bforest 来判断
    st_dbf, ps_dbf = setup(6, 8; Ta=20.0, θ0=0.3)
    st_crop, ps_crop = setup(25, 8; Ta=20.0, θ0=0.3)

    # 测试 ParamBEPS 结构存在
    @test ps_dbf.ψ_min isa Real
    @test ps_dbf.alpha isa Real
    @test ps_crop.ψ_min isa Real
    @test ps_crop.alpha isa Real
  end


  @testset "setup_jl" begin
    soil, st, ps = setup_jl(25, 8; Ta=20.0, θ0=0.3)

    @test soil isa Soil
    @test st isa StateBEPS
    @test ps isa ParamBEPS

    @test soil.n_layer == 5
    @test st.n_layer == 5
    @test ps.N == 5
  end


  @testset "setup_c" begin
    soil, state, ps = setup_c(25, 8; Ta=20.0, θ0=0.3)

    @test soil isa Soil_c
    @test state isa Vector{Float64}
    @test length(state) == 41
    @test ps isa ParamBEPS

    @test soil.n_layer == 5
  end

end
