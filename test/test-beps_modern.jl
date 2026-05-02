using BEPS, DataFrames, Test, Dates

LAI = readdlm(path_proj("examples/input/p1_LAI.txt"))[:]
forcing = deserialize(path_proj("data/p1_forcing"))
dates = DateTime(2010):Hour(1):DateTime(2010, 12, 31, 23)

VegType, SoilType = 25, 8
kw = (lon=120.5, lat=30.5, clumping=0.85, fix_snowpack=false)

ps = ParamBEPS(VegType, SoilType)
Ta = Float64(forcing.Tair[1])
state0, _ = setup(ps; Ta, Tsoil=2.2, θ0=0.4115, z_snow=0.0)

@testset "beps_modern" begin

  # --- 1. 基本运行 ---
  @testset "basic run" begin
    df, df_ET, states = beps_modern(forcing, LAI, dates; ps, state=state0, kw...)
    ntime = forcing.ntime

    @test size(df, 1) == ntime
    @test size(df_ET, 1) == ntime
    @test size(states.vectors.Tsoil_c) == (5, ntime)
    @test size(states.vectors.θ) == (5, ntime)
  end

  # --- 2. 物理边界 ---
  @testset "physical bounds" begin
    df, df_ET, states = beps_modern(forcing, LAI, dates; ps, state=state0, kw...)

    @test all(isfinite, states.vectors.Tsoil_c)
    @test all(isfinite, states.vectors.θ)
    @test all(states.vectors.θ .>= 0)
    @test all(states.vectors.θ .<= 1)
    @test all(isfinite, df.GPP)
    @test all(df.GPP .>= 0)
  end

  # --- 3. 参数变化影响输出 ---
  @testset "parameter sensitivity" begin
    df_ref, _, _ = beps_modern(forcing, LAI, dates; ps, state=state0, kw...)

    ps2 = deepcopy(ps)
    ps2.veg.VCmax25 *= 2.0
    df_vcmax, _, _ = beps_modern(forcing, LAI, dates; ps=ps2, state=state0, kw...)

    @test size(df_vcmax) == size(df_ref)
    @test sum(df_vcmax.GPP) > sum(df_ref.GPP)  # GPP 随 VCmax25 增大
  end

  # --- 4. state0 不被修改 ---
  @testset "state0 not mutated" begin
    θ_before = copy(state0.θ)
    beps_modern(forcing, LAI, dates; ps, state=state0, kw...)
    @test state0.θ == θ_before
  end

end
