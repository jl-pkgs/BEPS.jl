abstract type AbstractSoil end

@with_kw mutable struct Soil <: AbstractSoil
  flag        ::Cint    = Cint(0)
  n_layer     ::Cint    = Cint(5)
  step_period ::Cint    = Cint(1)

  Zp          ::Cdouble = Cdouble(0)
  Zsp         ::Cdouble = Cdouble(0)
  r_rain_g    ::Cdouble = Cdouble(0)
  soil_r      ::Cdouble = Cdouble(0)
  r_drainage  ::Cdouble = Cdouble(0)
  r_root_decay::Cdouble = Cdouble(0)
  ψ_min       ::Cdouble = Cdouble(0)
  alpha       ::Cdouble = Cdouble(0)
  f_soilwater ::Cdouble = Cdouble(0)

  dz          ::Vector{Float64} = zeros(10)
  f_root      ::Vector{Float64} = zeros(10)
  dt          ::Vector{Float64} = zeros(10)
  thermal_cond::Vector{Float64} = zeros(10)
  theta_vfc   ::Vector{Float64} = zeros(10)
  θ_vwp       ::Vector{Float64} = zeros(10)
  θ_sat       ::Vector{Float64} = zeros(10)
  Ksat        ::Vector{Float64} = zeros(10)
  ψ_sat       ::Vector{Float64} = zeros(10)
  b           ::Vector{Float64} = zeros(10) # Cambell parameter b
  density_soil::Vector{Float64} = zeros(10)
  f_org       ::Vector{Float64} = zeros(10) # organic matter
  ice_ratio   ::Vector{Float64} = zeros(10) # ice ratio
  θ           ::Vector{Float64} = zeros(10) # soil moisture
  θ_prev      ::Vector{Float64} = zeros(10) # soil moisture in previous time
  Tsoil_p     ::Vector{Float64} = zeros(10) # soil temperature in previous time
  Tsoil_c     ::Vector{Float64} = zeros(10) # soil temperature in current time

  f_water     ::Vector{Float64} = zeros(10) 
  ψ           ::Vector{Float64} = zeros(10)
  θb          ::Vector{Float64} = zeros(10) # not used, soil water content at the bottom of each layer
  ψb          ::Vector{Float64} = zeros(10) # not used
  r_waterflow ::Vector{Float64} = zeros(10)
  km          ::Vector{Float64} = zeros(10) # hydraulic conductivity
  Kb          ::Vector{Float64} = zeros(10) # not used
  KK          ::Vector{Float64} = zeros(10) # average conductivity of two soil layers
  Cs          ::Vector{Float64} = zeros(10)
  lambda      ::Vector{Float64} = zeros(10)
  Ett         ::Vector{Float64} = zeros(10)
  G           ::Vector{Float64} = zeros(10)

  # temporary variables in soil_water_factor_v2
  ft          ::Vector{Float64} = zeros(10)
  dtt         ::Vector{Float64} = zeros(10)
  fpsisr      ::Vector{Float64} = zeros(10)
end
