import Parameters: @with_kw, @with_kw_noshow

Value = getindex
Value! = setindex!

const FT = Cdouble

init_dbl() = Ref(0.0)
dbl() = Cdouble(0)



include("Constant.jl")
include("Leaf.jl")
include("Soil.jl")
include("Other_structs.jl")

include("InterTempVars.jl")
include("OutputET.jl")

export 
  Leaf, LeafRef,
  Soil, TSoil,
  ClimateData, Results, Cpools,
  InterTempVars,
  OutputET, Radiation

export FT, init_dbl, dbl
