module BEPS
# using BEPS

using UnPack
import Parameters: @with_kw, @with_kw_noshow
using Printf
using Reexport

export besp_main
export clang


# import Statistics: mean, std
include("IO.jl")
include("Ipaper.jl")
# include("c2julia.jl")

include("DataType/DataType.jl")
include("helpers.jl")
include("BEPS_helper.jl")

include("clang/BEPS_c.jl")
@reexport using BEPS.clang;

include("BEPS/BEPS.jl")

include("beps_inter_prg.jl")
include("beps_main.jl")

end # module BEPS
