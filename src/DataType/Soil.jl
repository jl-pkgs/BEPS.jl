abstract type AbstractSoil end

@with_kw mutable struct Soil <: AbstractSoil
  flag        ::Cint    = Cint(0)
  n_layer     ::Cint    = Cint(5)
  step_period ::Cint    = Cint(1)

  depth_water ::Cdouble = Cdouble(0)
  depth_snow  ::Cdouble = Cdouble(0)

  r_rain_g    ::Cdouble = Cdouble(0)
  soil_r      ::Cdouble = Cdouble(0)
  r_drainage  ::Cdouble = Cdouble(0)
  r_root_decay::Cdouble = Cdouble(0)
  ψ_min       ::Cdouble = Cdouble(0)
  alpha       ::Cdouble = Cdouble(0)
  f_soilwater ::Cdouble = Cdouble(0)

  dz          ::Vector{Float64} = zeros(10)
  f_root      ::Vector{Float64} = zeros(10) # root fraction
  dt          ::Vector{Float64} = zeros(10) # soil water stress factor
  thermal_cond::Vector{Float64} = zeros(10) # thermal conductivity
  theta_vfc   ::Vector{Float64} = zeros(10) # volumetric field capacity
  θ_vwp       ::Vector{Float64} = zeros(10) # volumetric wilting point
  θ_sat       ::Vector{Float64} = zeros(10) # volumetric saturation
  Ksat        ::Vector{Float64} = zeros(10) # saturated hydraulic conductivity
  ψ_sat       ::Vector{Float64} = zeros(10) # soil matric potential at saturation
  b           ::Vector{Float64} = zeros(10) # Cambell parameter b
  density_soil::Vector{Float64} = zeros(10) # soil density
  f_org       ::Vector{Float64} = zeros(10) # organic matter
  ice_ratio   ::Vector{Float64} = zeros(10) # ice ratio
  θ           ::Vector{Float64} = zeros(10) # soil moisture
  θ_prev      ::Vector{Float64} = zeros(10) # soil moisture in previous time
  Tsoil_p     ::Vector{Float64} = zeros(10) # soil temperature in previous time
  Tsoil_c     ::Vector{Float64} = zeros(10) # soil temperature in current time

  f_water     ::Vector{Float64} = zeros(10) # not used
  ψ           ::Vector{Float64} = zeros(10) # soil matric potential
  θb          ::Vector{Float64} = zeros(10) # not used, θ at the bottom of each layer
  ψb          ::Vector{Float64} = zeros(10) # not used
  r_waterflow ::Vector{Float64} = zeros(10) # vertical water flow rate
  km          ::Vector{Float64} = zeros(10) # hydraulic conductivity at middle point
  Kb          ::Vector{Float64} = zeros(10) # not used
  KK          ::Vector{Float64} = zeros(10) # average conductivity of two soil layers
  Cs          ::Vector{Float64} = zeros(10) 
  lambda      ::Vector{Float64} = zeros(10) 
  Ett         ::Vector{Float64} = zeros(10) # ET in each layer. derived var
  G           ::Vector{Float64} = zeros(10) 

  # temporary variables in soil_water_factor_v2
  ft          ::Vector{Float64} = zeros(10) 
  dtt         ::Vector{Float64} = zeros(10) 
  fpsisr      ::Vector{Float64} = zeros(10)
end
