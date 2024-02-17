@with_kw mutable struct Soil
  flag        ::Cint    = Cint(0)
  n_layer     ::Cint    = Cint(5)
  step_period ::Cint    = Cint(1)

  Zp          ::Cdouble = Cdouble(0)
  Zsp         ::Cdouble = Cdouble(0)
  r_rain_g    ::Cdouble = Cdouble(0)
  soil_r      ::Cdouble = Cdouble(0)
  r_drainage  ::Cdouble = Cdouble(0)
  r_root_decay::Cdouble = Cdouble(0)
  psi_min     ::Cdouble = Cdouble(0)
  alpha       ::Cdouble = Cdouble(0)
  f_soilwater ::Cdouble = Cdouble(0)

  d_soil      ::Vector{Float64} = zeros(10)
  f_root      ::Vector{Float64} = zeros(10)
  dt          ::Vector{Float64} = zeros(10)
  thermal_cond::Vector{Float64} = zeros(10)
  theta_vfc   ::Vector{Float64} = zeros(10)
  theta_vwp   ::Vector{Float64} = zeros(10)
  fei         ::Vector{Float64} = zeros(10)
  Ksat        ::Vector{Float64} = zeros(10)
  psi_sat     ::Vector{Float64} = zeros(10)
  b           ::Vector{Float64} = zeros(10) # Cambell parameter b
  density_soil::Vector{Float64} = zeros(10)
  f_org       ::Vector{Float64} = zeros(10) # organic matter
  ice_ratio   ::Vector{Float64} = zeros(10) # ice ratio
  thetam      ::Vector{Float64} = zeros(10) # soil moisture
  thetam_prev ::Vector{Float64} = zeros(10) # soil moisture in previous time
  temp_soil_p ::Vector{Float64} = zeros(10) # soil temperature in previous time
  temp_soil_c ::Vector{Float64} = zeros(10) # soil temperature in current time
  f_ice       ::Vector{Float64} = zeros(10) 
  psim        ::Vector{Float64} = zeros(10)
  thetab      ::Vector{Float64} = zeros(10) # soil water content at the bottom of each layer
  psib        ::Vector{Float64} = zeros(10)
  r_waterflow ::Vector{Float64} = zeros(10)
  km          ::Vector{Float64} = zeros(10) # hydraulic conductivity
  Kb          ::Vector{Float64} = zeros(10) 
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
