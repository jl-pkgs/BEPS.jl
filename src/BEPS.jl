module BEPS
# using BEPS

using UnPack
import Parameters: @with_kw, @with_kw_noshow
using Printf
using Reexport

export besp_main
export clang


path_proj(f...) = normpath(joinpath(@__DIR__, "..", f...))
libbeps = path_proj("deps/libbeps.dll")

# import Statistics: mean, std
include("IO.jl")
include("Ipaper.jl")
# include("c2julia.jl")

include("DataType/DataType.jl")
include("helpers.jl")
include("BEPS_helper.jl")

include("BEPS/BEPS.jl")

include("clang/BEPS_c.jl")
@reexport using BEPS.clang;

include("beps_inter_prg.jl")
include("beps_main.jl")

end # module BEPS
