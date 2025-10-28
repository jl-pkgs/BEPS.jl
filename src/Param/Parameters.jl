export VegParam
export ParamSoilHydraulic, ParamSoilThermal, ParamSoil,
  ParamSoilHydraulicLayers, ParamSoilThermalLayers
export BEPSmodel

using Parameters, DataFrames
import FieldMetadata: @metadata, @units, units
@metadata bounds nothing


@with_kw mutable struct VegParam{FT<:AbstractFloat}
  # lc::Int = 1
  # Ω0::FT = 0.7             # // clumping_index
  LAI_max_o::FT = 4.5     # ? 应该不重要
  LAI_max_u::FT = 2.4

  α_canopy_vis::FT = 0.055 # 0.04
  α_canopy_nir::FT = 0.300 # 0.25
  α_soil_sat::FT = 0.10    # albedo of saturated/dry soil, `rainfall1`
  α_soil_dry::FT = 0.35    # the albedo of dry soil

  # r_drainage::FT = 0.5     # ? 产流比例
  # r_root_decay::FT = 0.97  # ? decay_rate_of_root_distribution

  z_canopy_o::FT = 1.0      # [m]
  z_canopy_u::FT = 0.2      # [m]
  z_wind::FT = 2            # [m]

  g1_w::FT = 8             # Ball-Berry, slope coefficient
  g0_w::FT = 0.0175        # Ball-Berry, intercept_for_H2O

  VCmax25::FT = 89.45       # maximum capacity of Rubisco at 25℃
  # Jmax25::FT = 2.39 * 57.7 - 14.2 # 

  # # coefficient reflecting the sensitivity of stomata to VPD/moderately N-stressed plants
  N_leaf::FT = 1.74 + 0.71   # leaf Nitrogen content, mean value + 1 SD [g/m2]
  slope_Vc::FT = 33.79 / 57.7
end

include("utilize.jl")
include("readVegParam.jl")
include("macro.jl")


# 水力参数
@bounds @with_kw mutable struct ParamSoilHydraulic{FT<:AbstractFloat}
  # θ_vfc::FT = FT(0.3) | (0.1, 0.4)  # not used, volumetric field capacity
  θ_vwp::FT = FT(0.1) | (0.1, 0.5)  # volumetric wilting point
  θ_sat::FT = FT(0.45) | (0.1, 0.6) # volumetric saturation
  K_sat::FT = FT(1e-5) | (0.1, 0.7) # saturated hydraulic conductivity
  ψ_sat::FT = FT(-0.5) | (0.1, 0.9) # soil matric potential at saturation
  b::FT = FT(5.0) | (0.1, 0.4) # Cambell parameter b
end

# 热力参数
@with_kw mutable struct ParamSoilThermal{FT<:AbstractFloat}
  κ_dry::FT = FT(0.2)     # thermal conductivity
  ρ_soil::FT = FT(1300.0) # soil density, 
  V_SOM::FT = FT(0.02)    # organic matter fraction, Soil Organic Matter
end

@make_layers_struct ParamSoilHydraulic
@make_layers_struct ParamSoilThermal

@with_kw mutable struct ParamSoil{FT<:AbstractFloat}
  hydraulic::ParamSoilHydraulic{FT} = ParamSoilHydraulic{FT}()
  thermal::ParamSoilThermal{FT} = ParamSoilThermal{FT}()
end

@bounds @with_kw mutable struct BEPSmodel{FT<:AbstractFloat}
  N::Int = 5
  r_drainage::FT = Cdouble(0.50) | (0.2, 0.7)      # ? 地表排水速率（地表汇流），可考虑采用曼宁公式
  r_root_decay::FT = Cdouble(0.95) | (0.85, 0.999) # ? 根系分布衰减率, decay_rate_of_root_distribution

  ψ_min::FT = Cdouble(33.0)  # * 气孔关闭对应水势，33kPa，可根据植被类型指定
  alpha::FT = Cdouble(0.4)   # * 土壤水限制因子参数，He 2017 JGR-B, Eq. 4

  hydraulic::ParamSoilHydraulicLayers{FT} = ParamSoilHydraulicLayers{FT,N}()
  thermal::ParamSoilThermalLayers{FT} = ParamSoilThermalLayers{FT,N}()

  veg::VegParam{FT} = VegParam{FT}()
end


# 这里应该加一个show function，打印模型参数信息
function Base.show(io::IO, model::M) where {M<:BEPSmodel}
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
  println("-"^ss)
  printstyled(io, "Hydraulic: ", color=:blue, bold=true)
  print(io, model.hydraulic)

  println("-"^ss)
  printstyled(io, "Thermal: ", color=:blue, bold=true)
  print(io, model.thermal)

  println("-"^ss)
  printstyled(io, "Veg: ", color=:blue, bold=true)
  print(io, model.veg)
  print("-"^ss)
  return nothing
end


function init_soil!(soil::Soil, model::BEPSmodel{FT}) where {FT}
  # N = model.N
  soil.r_drainage = Cdouble(model.r_drainage)
  soil.r_root_decay = Cdouble(model.r_root_decay)
  soil.ψ_min = Cdouble(model.ψ_min)
  soil.alpha = Cdouble(model.alpha)

  soil.θ_vfc .= Cdouble(model.hydraulic.θ_vfc)
  soil.θ_vwp .= Cdouble(model.hydraulic.θ_vwp)
  soil.θ_sat .= Cdouble(model.hydraulic.θ_sat)
  soil.Ksat .= Cdouble(model.hydraulic.K_sat)
  soil.ψ_sat .= Cdouble(model.hydraulic.ψ_sat)
  soil.b .= Cdouble(model.hydraulic.b)

  soil.κ_dry .= Cdouble(model.thermal.κ_dry)
  soil.density_soil .= Cdouble(model.thermal.ρ_soil)
  soil.V_SOM .= Cdouble(model.thermal.V_SOM)
  return soil
end


# @with_kw mutable struct SoilParam{FT<:AbstractFloat}
#   n_layer::Cint = 5 # 土壤层数
#   # dz::Vector{Float64} = zeros(10) # 土壤厚度

#   r_drainage  ::FT = Cdouble(0.50)  # ? 地表排水速率（地表汇流），可考虑采用曼宁公式
#   r_root_decay::FT = Cdouble(0.95)  # ? 根系分布衰减率, decay_rate_of_root_distribution

#   ψ_min       ::FT = Cdouble(33.0)  # ? 气孔关闭对应水势，33kPa
#   alpha       ::FT = Cdouble(0.4)   # ? 土壤水限制因子参数，He 2017 JGR-B, Eq. 4

#   f_root::Vector{FT} = zeros(FT, 10) # * 根系比例，可根据r_root_decay设定

#   θ_vfc       ::Vector{FT} = zeros(FT, 10) # ? volumetric field capacity
#   θ_vwp       ::Vector{FT} = zeros(FT, 10) # ? volumetric wilting point
#   θ_sat       ::Vector{FT} = zeros(FT, 10) # ? volumetric saturation
#   Ksat        ::Vector{FT} = zeros(FT, 10) # ? saturated hydraulic conductivity
#   ψ_sat       ::Vector{FT} = zeros(FT, 10) # ? soil matric potential at saturation
#   b           ::Vector{FT} = zeros(FT, 10) # ? Cambell parameter b

#   κ_dry       ::Vector{FT} = zeros(FT, 10) # ? thermal conductivity
#   density_soil::Vector{FT} = zeros(FT, 10) # ? 土壤容重，soil density, for volume heat capacity
#   V_SOM       ::Vector{FT} = zeros(FT, 10) # ? 有机质含量，organic matter, for volume heat capacity
# end
