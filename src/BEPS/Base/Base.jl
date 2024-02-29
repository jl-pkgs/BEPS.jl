using DocStringExtensions: TYPEDFIELDS

include("Ipaper.jl")
include("helpers.jl")
include("s_coszs.jl")
include("lai2.jl")
include("aerodynamic_conductance.jl")

include("AirLayer.jl")

export s_coszs, lai2, aerodynamic_conductance_jl
