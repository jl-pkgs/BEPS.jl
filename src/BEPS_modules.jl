export s_coszs, lai2, ReadParamVeg,
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


# include("Base/Base.jl")
include("Param/Parameters.jl")

include("SoilPhysics/SoilPhysics.jl")
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


# 25: C3 农作物 (C3 Crops)
# 40: C4 农作物/草地 (C4 Crops/Grass)
function lai2!(veg::ParamVeg{T}, Ω::T, CosZs::T, lai::T,
  LAI::Leaf, PAI::Leaf) where {T<:AbstractFloat}

  lai_o = lai < 0.1 ? 0.1 : lai
  lai_u = !veg.has_understory ? 0.01 : 1.18 * exp(-0.99 * lai_o)
  lai_u > lai_o && (lai_u = 0.01)

  stem_o = veg.LAI_max_o * 0.2 # 
  stem_u = veg.LAI_max_u * 0.2
  lai2!(Ω, CosZs, stem_o, stem_u, lai_o, lai_u, LAI, PAI)
  lai_o, lai_u, stem_o, stem_u
end
