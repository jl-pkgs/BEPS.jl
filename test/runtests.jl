using Test, BEPS
using BEPS: path_proj

function nanmax(x)
  x = collect(x)
  x = x[.!isnan.(x)]
  maximum(x)
end

include("test-ModelParams.jl")
include("test-utilize.jl")
include("modules/modules.jl")
include("test-beps_main.jl")
