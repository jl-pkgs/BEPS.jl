import Parameters: @with_kw, @with_kw_noshow

Value = getindex
Value! = setindex!

const FT = Cdouble

# dbl() = Cdouble(0)
init_dbl() = Ref(0.0)
nzero(n) = tuple(zeros(n)...) # n double zero



include("Constant.jl")
include("Leaf.jl")
include("CanopyLayer.jl")
include("Soil.jl")

include("InterTempVars.jl")
include("Radiation.jl")
include("INPUT.jl")
include("OUTPUT.jl")


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
const TypeDF = Union{Results,ClimateData,OutputET}

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
  nothing
end


export 
  Leaf, 
  Soil, AbstractSoil, 
  ClimateData, Results, Cpools,
  InterTempVars,
  OutputET, Radiation

export FT, init_dbl
