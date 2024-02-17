module BEPS
# using BEPS

using UnPack
import Parameters: @with_kw, @with_kw_noshow
using Printf
using Reexport
@reexport using Serialization: deserialize, serialize
@reexport using DelimitedFiles: readdlm
export besp_main
export init_soil!
export Soil_c

path_proj(f...) = normpath(joinpath(@__DIR__, "..", f...))
libbeps = path_proj("deps/libbeps.dll")

# import Statistics: mean, std
include("DataFrames.jl")
include("Ipaper.jl")
# include("c2julia.jl")

include("DataType/DataType.jl")
include("helpers.jl")

include("BEPS/BEPS.jl")

include("clang/BEPS_c.jl")
@reexport import BEPS.clang;
import BEPS.clang: inter_prg_c, photosynthesis_c, Soil_c

include("BEPS_helper.jl")

include("beps_inter_prg.jl")
include("beps_main.jl")

end # module BEPS
