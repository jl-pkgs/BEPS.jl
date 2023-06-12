module BEPS
# using BEPS

using UnPack
import Parameters: @with_kw, @with_kw_noshow
using Printf

# import Statistics: mean, std
# using TimerOutputs
include("IO.jl")

# include("c2julia.jl")
include("Ipaper.jl")

include("clang/BEPS_c.jl")
# include("DataType.jl")
# include("BEPS_helper.jl")

# include("Soil.jl")
# include("init_soil.jl")

end # module BEPS
