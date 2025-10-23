abstract type AbstractSoil end

@with_kw mutable struct Soil <: AbstractSoil
  flag        ::Cint    = Cint(0)
  n_layer     ::Cint    = Cint(5)
  step_period ::Cint    = Cint(1)

  z_water ::Cdouble = Cdouble(0) # [state]
  z_snow  ::Cdouble = Cdouble(0) # [state]

  r_rain_g    ::Cdouble = Cdouble(0)        # the rainfall rate, un--on understory g--on ground surface  m/s
  soil_r      ::Cdouble = Cdouble(0)
  r_drainage  ::Cdouble = Cdouble(0)        # 土壤排水速率
  r_root_decay::Cdouble = Cdouble(0)        # decay_rate_of_root_distribution, 根系分布衰减率
  ψ_min       ::Cdouble = Cdouble(0)
  alpha       ::Cdouble = Cdouble(0)
  f_soilwater ::Cdouble = Cdouble(0)        # 土壤水限制因子

  dz          ::Vector{Float64} = zeros(10) # 土壤厚度
  f_root      ::Vector{Float64} = zeros(10) # 根系比例，root fraction
  dt          ::Vector{Float64} = zeros(10) # 土壤水限制因子，soil water stress factor
  κ_dry       ::Vector{Float64} = zeros(10) # thermal conductivity
  θ_vfc       ::Vector{Float64} = zeros(10) # volumetric field capacity
  θ_vwp       ::Vector{Float64} = zeros(10) # volumetric wilting point
  θ_sat       ::Vector{Float64} = zeros(10) # volumetric saturation
  Ksat        ::Vector{Float64} = zeros(10) # saturated hydraulic conductivity
  ψ_sat       ::Vector{Float64} = zeros(10) # soil matric potential at saturation
  b           ::Vector{Float64} = zeros(10) # Cambell parameter b
  density_soil::Vector{Float64} = zeros(10) # 土壤容重，soil density
  f_org       ::Vector{Float64} = zeros(10) # 有机质含量，organic matter
  ice_ratio   ::Vector{Float64} = zeros(10) # [state]，ice ratio，
  θ           ::Vector{Float64} = zeros(10) # [state], soil moisture
  θ_prev      ::Vector{Float64} = zeros(10) # [state], soil moisture in previous time
  Tsoil_p     ::Vector{Float64} = zeros(10) # [state], soil temperature in previous time
  Tsoil_c     ::Vector{Float64} = zeros(10) # [state], soil temperature in current time

  f_water     ::Vector{Float64} = zeros(10) # not used
  ψ           ::Vector{Float64} = zeros(10) # [state], soil matric potential
  θb          ::Vector{Float64} = zeros(10) # not used, θ at the bottom of each layer
  ψb          ::Vector{Float64} = zeros(10) # not used
  r_waterflow ::Vector{Float64} = zeros(10) # [state], vertical water flow rate
  km          ::Vector{Float64} = zeros(10) # [state], hydraulic conductivity at middle point
  Kb          ::Vector{Float64} = zeros(10) # not used
  KK          ::Vector{Float64} = zeros(10) # [state], average conductivity of two soil layers
  Cs          ::Vector{Float64} = zeros(10) # [state], volume heat capacity
  κ           ::Vector{Float64} = zeros(10) # [state]
  Ett         ::Vector{Float64} = zeros(10) # 每层蒸发量ET in each layer
  G           ::Vector{Float64} = zeros(10) # 土壤热通量

              # temporary variables in soil_water_factor_v2
  ft          ::Vector{Float64} = zeros(10) 
  dtt         ::Vector{Float64} = zeros(10) 
  fpsisr      ::Vector{Float64} = zeros(10)
end
