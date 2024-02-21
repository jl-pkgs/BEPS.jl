nzero(n) = tuple(zeros(n)...) # n double zero
const NT10 = NTuple{10,Cdouble}

@with_kw mutable struct Soil_c
  flag        ::Cint = Cint(0)
  n_layer     ::Cint = Cint(5)
  step_period ::Cint = Cint(1)
  Zp          ::Cdouble = Cdouble(0)
  Zsp         ::Cdouble = Cdouble(0)
  r_rain_g    ::Cdouble = Cdouble(0)
  soil_r      ::Cdouble = Cdouble(0)
  r_drainage  ::Cdouble = Cdouble(0)
  r_root_decay::Cdouble = Cdouble(0)
  ψ_min     ::Cdouble = Cdouble(0)
  alpha       ::Cdouble = Cdouble(0)
  f_soilwater ::Cdouble = Cdouble(0)
  dz      ::NT10 = nzero(10)
  f_root      ::NT10 = nzero(10)
  dt          ::NT10 = nzero(10)
  thermal_cond::NT10 = nzero(10)
  theta_vfc   ::NT10 = nzero(10)
  theta_vwp   ::NT10 = nzero(10)
  θ_sat         ::NT10 = nzero(10)
  Ksat        ::NT10 = nzero(10)
  ψ_sat     ::NT10 = nzero(10)
  b           ::NT10 = nzero(10)
  density_soil::NT10 = nzero(10)
  f_org       ::NT10 = nzero(10)
  ice_ratio   ::NT10 = nzero(10)
  θ      ::NT10 = nzero(10)
  θ_prev ::NT10 = nzero(10)
  Tsoil_p ::NT10 = nzero(10)
  Tsoil_c ::NT10 = nzero(10)
  f_water       ::NT10 = nzero(10)
  ψ        ::NT10 = nzero(10)
  θb      ::NT10 = nzero(10)
  ψb        ::NT10 = nzero(10)
  r_waterflow ::NT10 = nzero(10)
  km          ::NT10 = nzero(10)
  Kb          ::NT10 = nzero(10)
  KK          ::NT10 = nzero(10)
  Cs          ::NT10 = nzero(10)
  lambda      ::NT10 = nzero(10)
  Ett         ::NT10 = nzero(10)
  G           ::NT10 = nzero(10)
end
