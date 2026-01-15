export ParamVeg
export ParamSoilHydraulic, ParamSoilThermal, ParamSoil,
  ParamSoilHydraulicLayers, ParamSoilThermalLayers
export ParamBEPS


using Parameters, DataFrames
import FieldMetadata: @metadata, @units, units
@metadata bounds nothing


include("GlobalData.jl")
include("macro.jl")

@bounds @with_kw mutable struct ParamVeg{FT<:AbstractFloat}
  # lc::Int = 1
  # Ω0::FT = 0.7             # // clumping_index
  has_understory::Bool = true      # 
  is_bforest::Bool = false         # broadleaf forest

  LAI_max_o::FT = 4.5 | (0.1, 7.0)      # LAI max for overstory
  LAI_max_u::FT = 2.4 | (0.1, 7.0)      # LAI max for understory

  α_canopy_vis::FT = 0.055 | (0.02, 0.15)  # canopy albedo visible
  α_canopy_nir::FT = 0.300 | (0.15, 0.50)  # canopy albedo near-infrared
  α_soil_sat::FT = 0.10 | (0.05, 0.20)     # albedo of saturated soil
  α_soil_dry::FT = 0.35 | (0.20, 0.50)     # albedo of dry soil

  # r_drainage::FT = 0.5     # ? 产流比例
  r_root_decay::FT = Cdouble(0.95) | (0.85, 0.999) # ? 根系分布衰减率, decay_rate_of_root_distribution

  z_canopy_o::FT = 1.0 | (0.1, 50.0)    # overstory canopy height [m]
  z_canopy_u::FT = 0.2 | (0.05, 5.0)    # understory canopy height [m]
  z_wind::FT = 2 | (1.0, 100.0)         # wind measurement height [m]

  g1_w::FT = 8 | (1.0, 20.0)            # Ball-Berry slope coefficient
  g0_w::FT = 0.0175 | (0.001, 0.1)      # Ball-Berry intercept for H2O

  VCmax25::FT = 89.45 | (5.0, 200.0)    # max Rubisco capacity at 25℃ [μmol m-2 s-1], global range
  # Jmax25::FT = 2.39 * 57.7 - 14.2 #

  # # coefficient reflecting the sensitivity of stomata to VPD/moderately N-stressed plants
  N_leaf::FT = 1.74 + 0.71 | (0.5, 5.0)   # leaf Nitrogen content, mean value + 1 SD [g/m2]
  slope_Vc::FT = 33.79 / 57.7 | (0.3, 1.0) # slope for Vcmax-N relationship
end


# 水力参数
@bounds @with_kw mutable struct ParamSoilHydraulic{FT<:AbstractFloat}
  θ_vfc::FT = FT(0.3) | (0.1, 0.4)  # volumetric field capacity
  θ_vwp::FT = FT(0.1) | (0.1, 0.5)  # volumetric wilting point
  θ_sat::FT = FT(0.45) | (0.1, 0.6) # volumetric saturation
  K_sat::FT = FT(1e-5) | (0.1, 0.7) # saturated hydraulic conductivity
  ψ_sat::FT = FT(-0.5) | (0.1, 0.9) # soil matric potential at saturation
  b::FT = FT(5.0) | (0.1, 0.4) # Cambell parameter b
end

# 热力参数
@bounds @with_kw mutable struct ParamSoilThermal{FT<:AbstractFloat}
  κ_dry::FT = FT(0.2) | (0.05, 0.5)          # dry soil thermal conductivity [W m-1 K-1]
  ρ_soil::FT = FT(1300.0) | (1000.0, 2000.0) # soil bulk density [kg m-3]
  V_SOM::FT = FT(0.02) | (0.0, 0.3)          # organic matter volume fraction [-]
end

@make_layers_struct ParamSoilHydraulic
@make_layers_struct ParamSoilThermal

@with_kw mutable struct ParamSoil{FT<:AbstractFloat}
  hydraulic::ParamSoilHydraulic{FT} = ParamSoilHydraulic{FT}()
  thermal::ParamSoilThermal{FT} = ParamSoilThermal{FT}()
end


@bounds @with_kw_noshow mutable struct ParamBEPS{FT<:AbstractFloat}
  N::Int = 5
  dz::Vector{FT} = FT[0.05, 0.10, 0.20, 0.40, 1.25]  # 土壤层厚度 [m], BEPS V2023
  r_drainage::FT = Cdouble(0.50) | (0.2, 0.7)  # ? 地表排水速率（地表汇流），可考虑采用曼宁公式

  ψ_min::FT = Cdouble(33.0)  # [m], about 0.10~0.33 MPa开始胁迫点
  alpha::FT = Cdouble(0.4)   # [-], 土壤水限制因子参数，He 2017 JGR-B, Eq. 4

  hydraulic::ParamSoilHydraulicLayers{FT} = ParamSoilHydraulicLayers{FT,N}()
  thermal::ParamSoilThermalLayers{FT} = ParamSoilThermalLayers{FT,N}()

  veg::ParamVeg{FT} = ParamVeg{FT}()
end

# `kw...`: other params like, `r_drainage`
function ParamBEPS(VegType::Int, SoilType::Int; N::Int=5, FT=Float64, kw...)
  veg = InitParam_Veg(VegType; FT)  
  hydraulic, thermal = InitParam_Soil(SoilType, N, FT)

  ψ_min = veg.is_bforest ? FT(10.0) : FT(33.0) # 开始胁迫点
  alpha = veg.is_bforest ? FT(1.5) : FT(0.4)   # 土壤水限制因子参数，He 2017 JGR-B, Eq. 4

  ParamBEPS{FT}(;
    N,  kw..., ψ_min, alpha,
    hydraulic, thermal, veg
  )
end


# 这里应该加一个show function，打印模型参数信息
function Base.show(io::IO, model::M) where {M<:ParamBEPS}
  printstyled(io, "$M, N = $(model.N)\n", color=:blue, bold=true)

  fields_all = fieldnames(M)
  fields = setdiff(fields_all, [:N, :hydraulic, :thermal, :veg])

  n = length(fields)
  for i = 1:n
    field = fields[i]
    value = getfield(model, field)
    type = typeof(value)
    isa(value, Function) && (type = Function)
    println(io, "  $field\t: {$type} $value")
    # (i != n) && print(io, "\n")
  end

  ss = 60
  println(io, "-"^ss)
  printstyled(io, "Hydraulic: ", color=:blue, bold=true)
  print(io, model.hydraulic)

  println(io, "-"^ss)
  printstyled(io, "Thermal: ", color=:blue, bold=true)
  print(io, model.thermal)

  println(io, "-"^ss)
  printstyled(io, "Veg: ", color=:blue, bold=true)
  print(io, model.veg)
  print("-"^ss)
  return nothing
end


include("InitParam.jl")
include("Sync.jl")
include("deprecated/VegHelper.jl")
include("deprecated/ReadParamVeg.jl")
include("deprecated/Init_Soil_Parameters.jl")
