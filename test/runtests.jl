using Test, BEPS
using BEPS: path_proj

function nanmax(x)
  x = collect(x)
  x = x[.!isnan.(x)]
  maximum(x)
end

include("test-photosynthesis_standalone.jl")
include("test-beps_main.jl")
include("test-beps_modern.jl")
include("test-ModelParams.jl")
include("test-utilize.jl")
include("test-optimization.jl")
include("test-evaluate.jl")
include("test-multisite.jl")
include("test-site-io.jl")
include("modules/modules.jl")
