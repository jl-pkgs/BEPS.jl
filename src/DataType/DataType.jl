import Parameters: @with_kw, @with_kw_noshow
const FT = Cdouble

Value = getindex
Value! = setindex!

# dbl() = Cdouble(0)
init_dbl() = Ref(0.0)
init_dbl(x::T) where {T<:Real} = Ref(x)
nzero(n) = tuple(zeros(n)...) # n double zero


include("Constant.jl")
# include("Leaf.jl")
include("CanopyLayer.jl")
include("BEPS_State.jl")
include("Params/Params.jl")

include("Met.jl")
include("macro.jl")
include("OUTPUT.jl")
include("setup.jl")


@with_kw mutable struct Radiation <: AbstractFlux
  Rs_o_df::FT = 0.0
  Rs_u_df::FT = 0.0

  Rs_o_dir::FT = 0.0
  Rs_u_dir::FT = 0.0

  Rns_o_df::FT = 0.0
  Rns_u_df::FT = 0.0
  Rns_g_df::FT = 0.0

  Rns_o_dir::FT = 0.0
  Rns_u_dir::FT = 0.0
  Rns_g_dir::FT = 0.0

  Rs_df::FT = 0.0
  Rs_dir::FT = 0.0
end

@with_kw mutable struct Cpools
  Ccd::NTuple{3,Cdouble} = nzero(3)
  Cssd::NTuple{3,Cdouble} = nzero(3)
  Csmd::NTuple{3,Cdouble} = nzero(3)
  Cfsd::NTuple{3,Cdouble} = nzero(3)
  Cfmd::NTuple{3,Cdouble} = nzero(3)
  Csm::NTuple{3,Cdouble} = nzero(3)
  Cm::NTuple{3,Cdouble} = nzero(3)
  Cs::NTuple{3,Cdouble} = nzero(3)
  Cp::NTuple{3,Cdouble} = nzero(3)
end


# # current not used
# @with_kw mutable struct TSoil
#   T_ground::Cdouble = 0.0
#   T_any0::Cdouble = 0.0
#   T_soil0::Cdouble = 0.0
#   T_snow::Cdouble = 0.0
#   T_snow1::Cdouble = 0.0
#   T_snow2::Cdouble = 0.0
#   G::Cdouble = 0.0
# end
# export TSoil


## fill valuesFlux
const TypeDF = Union{Flux,Met,ETFlux}

## put struct into a data.frame
function Base.getindex(x::T, i::Int)::FT where {T<:TypeDF}
  # key = fieldnames(T)[i]
  getfield(x, i)
end

Base.length(x::T) where {T<:TypeDF} = fieldcount(T)

@generated function fill_res!(df::DataFrame, res::T, k::Int) where {T<:TypeDF}
  fs = fieldnames(T)
  assigns = Vector{Any}(undef, length(fs))
  for (i, f) in pairs(fs)
    # QuoteNode 让字段名作为常量符号 → DataFrames 列直接定位 + 静态 getfield
    assigns[i] = :(@inbounds df[!, $(QuoteNode(f))][k] = getfield(res, $(QuoteNode(f))))
  end
  return Expr(:block, assigns..., :(return nothing))
end


export Leaf, Soil, AbstractSoil, 
  Met, Flux, Cpools, ETFlux, Radiation

export FT, init_dbl, set!
