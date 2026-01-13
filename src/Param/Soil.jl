abstract type AbstractSoil end

# TODO: 这里应该把状态变量与模型参数分隔开
# state, params = setup(model)

# ?     : 需要优化的参数
# state : 状态变量
# //    : 未使用的参数
@with_kw mutable struct Soil <: AbstractSoil
  flag        ::Cint    = Cint(0) # // not used
  n_layer     ::Cint    = Cint(5) # 土壤层数
  step_period ::Cint    = Cint(1) # // not used

  z_water ::Cdouble = Cdouble(0) # [state]
  z_snow  ::Cdouble = Cdouble(0) # [state]

  # the rainfall rate, un--on understory on ground surface  m/s
  r_rain_g    ::Cdouble = Cdouble(0)        # [state], 达到地地表降水, PE, [m/s]

  soil_r      ::Cdouble = Cdouble(0)        # // not used, soil surface resistance for water
  r_drainage  ::Cdouble = Cdouble(0)        # ? 地表排水速率（地表汇流）
  r_root_decay::Cdouble = Cdouble(0)        # ? 根系分布衰减率, decay_rate_of_root_distribution
  ψ_min       ::Cdouble = Cdouble(0)        # ? 开始胁迫，33[m] = 0.33[MPa]
  alpha       ::Cdouble = Cdouble(0)        # ? 土壤水限制因子参数，He 2017 JGR-B, Eq. 4
  f_soilwater ::Cdouble = Cdouble(0)        # [state], 总体的土壤水限制因子

  dz          ::Vector{Float64} = zeros(10) # 土壤厚度
  f_root      ::Vector{Float64} = zeros(10) # [state], 根系比例，root fraction
  dt          ::Vector{Float64} = zeros(10) # [state], 每层的土壤水限制因子，已归一化
  κ_dry       ::Vector{Float64} = zeros(10) # ? thermal conductivity
  θ_vfc       ::Vector{Float64} = zeros(10) # ? volumetric field capacity
  θ_vwp       ::Vector{Float64} = zeros(10) # ? volumetric wilting point
  θ_sat       ::Vector{Float64} = zeros(10) # ? volumetric saturation
  Ksat        ::Vector{Float64} = zeros(10) # ? saturated hydraulic conductivity
  ψ_sat       ::Vector{Float64} = zeros(10) # ? soil matric potential at saturation
  b           ::Vector{Float64} = zeros(10) # ? Cambell parameter b
  ρ_soil      ::Vector{Float64} = zeros(10) # ? 土壤容重，soil density, for volume heat capacity
  V_SOM       ::Vector{Float64} = zeros(10) # ? 有机质含量，organic matter, for volume heat capacity

  ice_ratio   ::Vector{Float64} = zeros(10) # [state]，ice ratio，
  θ           ::Vector{Float64} = zeros(10) # [state], soil moisture
  θ_prev      ::Vector{Float64} = zeros(10) # [state], soil moisture in previous time
  Tsoil_p     ::Vector{Float64} = zeros(10) # [state], soil temperature in previous time
  Tsoil_c     ::Vector{Float64} = zeros(10) # [state], soil temperature in current time

  f_water     ::Vector{Float64} = zeros(10) # // not used
  ψ           ::Vector{Float64} = zeros(10) # [state], soil matric potential
  θb          ::Vector{Float64} = zeros(10) # // not used, θ at the bottom of each layer
  ψb          ::Vector{Float64} = zeros(10) # // not used
  r_waterflow ::Vector{Float64} = zeros(10) # [state], vertical water flow rate
  km          ::Vector{Float64} = zeros(10) # [state], hydraulic conductivity at middle point
  Kb          ::Vector{Float64} = zeros(10) # // not used
  KK          ::Vector{Float64} = zeros(10) # [state], average conductivity of two soil layers
  Cs          ::Vector{Float64} = zeros(10) # [state], volume heat capacity
  κ           ::Vector{Float64} = zeros(10) # [state]
  Ett         ::Vector{Float64} = zeros(10) # [state], 每层蒸发量ET in each layer
  G           ::Vector{Float64} = zeros(10) # [state], 土壤热通量

  ## temporary variables in soil_water_factor_v2
  ft          ::Vector{Float64} = zeros(10) # [state], f_i(Tsoil_i), 温度对水分限制影响, Eq. 5
  dtt         ::Vector{Float64} = zeros(10) # [state], 叠加根系分布比例，f_root[i] * fpsisr[i]
  fpsisr      ::Vector{Float64} = zeros(10) # [state], f_{w,i}, He et al., 2017, Eq. 3
end


@with_kw mutable struct State{FT}
  "Surface Temperature: [T_ground, T_surf_snow, T_surf_mix, T_snow_L1, T_snow_L2]"
  Ts::Vector{FT} = zeros(FT, 5)         # 4:8
  Ts_prev::Vector{FT} = zeros(FT, 5) # 10:15
  θ_prev::Vector{FT} = zeros(FT, 5)      # 22:27
  ice_ratio::Vector{FT} = zeros(FT, 5)   # 28:33

  Qhc_o::FT = 0.0                    # [11] sensible heat flux
  m_water::Layer2 = Layer2{FT}()     # [15, 18] + 1
  m_snow::Layer3 = Layer3{FT}()      # [16, 19, 20] + 1
  ρ_snow::FT = 250.0                 # [kg m-3] snow density
end
# 拖着`ρ_snow`，`ρ_snow`也是一个状态连续的变量
# https://www.eoas.ubc.ca/courses/atsc113/snow/met_concepts/07-met_concepts/07b-newly-fallen-snow-density/

export State, Soil
