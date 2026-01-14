#import "@local/modern-cug-report:0.1.3": *
#show: doc => template(doc, footer: "CUG水文气象学2025", header: "")

#set page(margin: (
  top: 1cm,
  // bottom: 2cm,
  x: 1.5cm,
))

*BEPS模型*：*参数统一架构*
#v(-0.6em)

```julia
@with_kw mutable struct VegParam{FT<:AbstractFloat}
  # lc::Int = 1
  # Ω0::FT = 0.7            # // clumping_index
  LAI_max_o::FT = 4.5       # ? 应该不重要
  LAI_max_u::FT = 2.4

  α_canopy_vis::FT = 0.055  # 0.04
  α_canopy_nir::FT = 0.300  # 0.25
  α_soil_sat::FT = 0.10     # albedo of saturated/dry soil, `rainfall1`
  α_soil_dry::FT = 0.35     # the albedo of dry soil

  # r_drainage::FT = 0.5     # ? 产流比例
  # r_root_decay::FT = 0.97  # ? decay_rate_of_root_distribution
  z_canopy_o::FT = 1.0      # [m]
  z_canopy_u::FT = 0.2      # [m]
  z_wind::FT = 2            # [m]

  g1_w::FT = 8              # Ball-Berry, slope coefficient
  g0_w::FT = 0.0175         # Ball-Berry, intercept_for_H2O

  VCmax25::FT = 89.45       # maximum capacity of Rubisco at 25℃
  # Jmax25::FT = 2.39 * 57.7 - 14.2 # 

  ## coefficient reflecting the sensitivity of stomata to VPD/moderately N-stressed plants
  N_leaf::FT = 1.74 + 0.71  # leaf Nitrogen content, mean value + 1 SD [g/m2]
  slope_Vc::FT = 33.79 / 57.7
end

@with_kw mutable struct SoilParam{FT<:AbstractFloat}
  n_layer::Cint = 5 # 土壤层数
  # dz::Vector{Float64} = zeros(10) # 土壤厚度

  r_drainage  ::FT = Cdouble(0.50)  # ? 地表排水速率（地表汇流），可考虑采用曼宁公式
  r_root_decay::FT = Cdouble(0.95)  # ? 根系分布衰减率, decay_rate_of_root_distribution

  ψ_min       ::FT = Cdouble(33.0)  # ? 气孔关闭对应水势，33kPa
  alpha       ::FT = Cdouble(0.4)   # ? 土壤水限制因子参数，He 2017 JGR-B, Eq. 4

  f_root::Vector{FT} = zeros(FT, 10) # * 根系比例，可根据r_root_decay设定

  θ_vfc       ::Vector{FT} = zeros(FT, 10) # ? volumetric field capacity
  θ_vwp       ::Vector{FT} = zeros(FT, 10) # ? volumetric wilting point
  θ_sat       ::Vector{FT} = zeros(FT, 10) # ? volumetric saturation
  Ksat        ::Vector{FT} = zeros(FT, 10) # ? saturated hydraulic conductivity
  ψ_sat       ::Vector{FT} = zeros(FT, 10) # ? soil matric potential at saturation
  b           ::Vector{FT} = zeros(FT, 10) # ? Cambell parameter b

  ## for volume heat capacity
  κ_dry       ::Vector{FT} = zeros(FT, 10) # ? thermal conductivity
  density_soil::Vector{FT} = zeros(FT, 10) # ? 土壤容重，soil density
  f_org       ::Vector{FT} = zeros(FT, 10) # ? 有机质含量，organic matter
end
```
