import Parameters: @with_kw, @with_kw_noshow
const FT = Cdouble

Value = getindex
Value! = setindex!

# dbl() = Cdouble(0)
init_dbl() = Ref(0.0)
nzero(n) = tuple(zeros(n)...) # n double zero


include("Constant.jl")
include("Leaf.jl")
include("CanopyLayer.jl")
include("Soil.jl")

include("InterTempVars.jl")
include("Met.jl")
include("OUTPUT.jl")


@with_kw mutable struct Radiation
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

@with_kw mutable struct State{FT}
  "Surface Temperature: [Ts0, Tsn, Tsm0, Tsn1, Tsn2]"
  Ts::Vector{FT} = zeros(FT, 5)         # 4:8
  Tsoil_prev::Vector{FT} = zeros(FT, 6) # 10:15
  θ_prev::Vector{FT} = zeros(FT, 6)      # 22:27
  ice_ratio::Vector{FT} = zeros(FT, 6)   # 28:33

  Qhc_o::FT = 0.0
  m_water::Layer2 = Layer2{FT}()     # [15, 18] + 1
  m_snow::Layer3 = Layer3{FT}() # [16, 19, 20] + 1
end

export State

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


## fill values
const TypeDF = Union{Results,Met,OutputET}

## put struct into a data.frame
function Base.getindex(x::T, i::Int)::FT where {T<:TypeDF}
  # key = fieldnames(T)[i]
  getfield(x, i)
end

Base.length(x::T) where {T<:TypeDF} = fieldcount(T)

function fill_res!(df::DataFrame, Res::T, k::Int) where {T<:TypeDF}
  n = length(Res)
  for i in 1:n
    df[k, i] = Res[i]
  end
  return nothing
end


export Leaf, Soil, AbstractSoil, 
  Met, Results, Cpools, InterTempVars, OutputET, Radiation

export FT, init_dbl, set!
