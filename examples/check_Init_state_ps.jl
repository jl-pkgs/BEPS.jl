using BEPS
using Test

const INIT_KW = (
  VegType=25,
  SoilType=8,
  Ta=10.0,
  Tsoil=2.2,
  θ0=0.4115,
  z_snow=0.0,
  r_drainage=0.5,
  r_root_decay=0.95,
)

is_scalar(x) = x isa Number || x isa AbstractString || x isa Symbol || x === nothing
is_same(a, b) = typeof(a) == typeof(b) && a == b

function push_diff!(diffs, path, a, b)
  if typeof(a) != typeof(b)
    push!(diffs, (; path, main_type=typeof(a), modern_type=typeof(b), main=a, modern=b))
    return
  end

  if is_scalar(a)
    is_same(a, b) || push!(diffs, (; path, main=a, modern=b))
    return
  end

  if a isa AbstractArray
    size(a) == size(b) || return push!(diffs, (; path, main_size=size(a), modern_size=size(b)))

    for i in eachindex(a, b)
      is_same(a[i], b[i]) || push!(diffs, (; path="$(path)[$i]", main=a[i], modern=b[i]))
    end
    return
  end

  for name in fieldnames(typeof(a))
    push_diff!(diffs, "$path.$name", getfield(a, name), getfield(b, name))
  end
end

function diff(a, b, root)
  diffs = NamedTuple[]
  push_diff!(diffs, root, a, b)
  diffs
end

init_kw(; kw...) = (; INIT_KW..., kw...)

function build_main_init(; kw...)
  init = init_kw(; kw...)
  (; VegType, SoilType, Ta, Tsoil, θ0, z_snow, r_drainage, r_root_decay) = init

  _, state, ps = setup_model(VegType, SoilType;
    version="julia", Ta, Tsoil, θ0, z_snow, r_drainage, r_root_decay)

  state, ps
end

function build_modern_init(; kw...)
  init = init_kw(; kw...)
  (; VegType, SoilType, Ta, Tsoil, θ0, z_snow, r_drainage, r_root_decay) = init

  ps = ParamBEPS(VegType, SoilType; r_drainage)
  ps.veg.r_root_decay = r_root_decay
  state, ps = setup(ps; Ta, Tsoil, θ0, z_snow)

  state, ps
end

function print_diffs(title, diffs)
  println("\n$title")
  isempty(diffs) && return println("  OK")

  for d in diffs
    println("  ", d)
  end
end

state_main, ps_main = build_main_init()
state_modern, ps_modern = build_modern_init()

state_diffs = diff(state_main, state_modern, "state")
ps_diffs = diff(ps_main, ps_modern, "ps")

print_diffs("StateBEPS init differences", state_diffs)
print_diffs("ParamBEPS init differences", ps_diffs)

@testset "beps_main and beps_modern init" begin
  @test isempty(state_diffs)
  @test isempty(ps_diffs)
end
