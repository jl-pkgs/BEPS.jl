import Parameters: @with_kw, @with_kw_noshow

Value = getindex
Value! = setindex!

const FT = Cdouble

# dbl() = Cdouble(0)
init_dbl() = Ref(0.0)
nzero(n) = tuple(zeros(n)...) # n double zero



include("Constant.jl")
include("Leaf.jl")
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


export 
  Leaf, 
  Soil, 
  ClimateData, Results, Cpools,
  InterTempVars,
  OutputET, Radiation

export FT, init_dbl
