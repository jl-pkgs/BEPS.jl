using BEPS, Test
import Parameters: @with_kw

# ── 测试用的最小 Flux / State 结构 ────────────────────────────────────────────
@with_kw mutable struct MockFlux <: AbstractFlux
  a::Float64 = 0.0
  b::Float64 = 0.0
end

@DefState MockState (:x, :y)
@DefFluxSeries MockSeries = MockFlux
@DefFluxSeries ExtraSeries = MockFlux extra_c
@DefStateSeries MockStateSer = MockState

@testset "DefFluxSeries / DefStateSeries macros" begin

  # ── 1. 生成结构的字段和类型 ──────────────────────────────────────────────────
  @testset "struct fields" begin
    s = MockSeries{Float64}(; ntime=10)
    @test s.ntime == 10
    @test s.a isa Vector{Float64}
    @test s.b isa Vector{Float64}
    @test length(s.a) == 10

    ss = MockStateSer{Float64}(; ntime=5)
    @test ss.ntime == 5
    @test ss.x isa Vector{Float64}
  end

  # ── 2. setindex!：series[i] = struct ─────────────────────────────────────────
  @testset "setindex!" begin
    s = MockSeries{Float64}(; ntime=3)
    f = MockFlux(a=1.5, b=2.5)
    s[2] = f
    @test s.a[2] == 1.5
    @test s.b[2] == 2.5
    @test s.a[1] == 0.0  # 未赋值的位置保持零
  end

  # ── 3. getindex 切片：series[r] ───────────────────────────────────────────────
  @testset "getindex slice" begin
    s = MockSeries{Float64}(; ntime=10)
    s.a .= 1:10
    s.b .= 11:20

    sub = s[3:6]
    @test sub.ntime == 4
    @test sub.a == [3.0, 4.0, 5.0, 6.0]
    @test sub.b == [13.0, 14.0, 15.0, 16.0]
  end

  # ── 4. AbstractStateSeries 的 getindex[i] 重建 State ─────────────────────────
  @testset "StateSeries getindex reconstruction" begin
    ss = MockStateSer{Float64}(; ntime=5)
    ss.x .= 1:5
    ss.y .= 6:10

    st = ss[3]
    @test st isa MockState{Float64}
    @test st.x == 3.0
    @test st.y == 8.0
  end

  # ── 5. 额外字段 ───────────────────────────────────────────────────────────────
  @testset "extra fields" begin
    s = ExtraSeries{Float64}(; ntime=4)
    @test hasproperty(s, :extra_c)
    @test length(s.extra_c) == 4
  end

  # ── 6. 超类型正确 ─────────────────────────────────────────────────────────────
  @testset "supertypes" begin
    @test MockSeries{Float64} <: AbstractFluxSeries{Float64}
    @test MockStateSer{Float64} <: AbstractStateSeries{Float64}
  end

end
