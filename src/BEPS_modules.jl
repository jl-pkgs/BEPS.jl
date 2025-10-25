# include("Base/Base.jl")
include("Param/Parameters.jl")

include("Soil/Soil.jl")
# include("soil_thermal_regime.jl")

include("aerodynamic_conductance.jl")
# include("latent_heat.jl")
# include("sensible_heat.jl")
include("heat_H_and_LE.jl")

# include("Leaf_Temperature.jl")
include("surface_temperature.jl")
# include("transpiration.jl")
include("evaporation_soil.jl")
include("evaporation_canopy.jl")

include("rainfall_stage.jl")
include("snowpack.jl")

include("netRadiation.jl")
include("photosynthesis.jl")

include("inter_prg.jl")

export s_coszs, lai2, readVegParam, 
  rainfall_stage1_jl, rainfall_stage2_jl, 
  snowpack_stage1_jl, snowpack_stage2_jl, snowpack_stage3_jl

export aerodynamic_conductance_jl, 
  sensible_heat_jl,
  latent_heat!,
  Leaf_Temperatures_jl, 
  surface_temperature_jl,
  transpiration_jl, 
  evaporation_canopy_jl,
  evaporation_soil_jl, 
  netRadiation_jl, 
  photosynthesis_jl
