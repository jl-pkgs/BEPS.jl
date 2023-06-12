import Parameters: @with_kw, @with_kw_noshow

Value = getindex
Value! = setindex!

init_dbl() = Ref(0.0)
dbl() = Cdouble(0)

# n double zero
nzero(n) = tuple(zeros(n)...)
const NT10 = NTuple{10,Cdouble}

@with_kw mutable struct Soil
  flag::Cint = Cint(0)
  n_layer::Cint = Cint(5)
  step_period::Cint = Cint(1)
  Zp::Cdouble = Cdouble(0)
  Zsp::Cdouble = Cdouble(0)
  r_rain_g::Cdouble = Cdouble(0)
  soil_r::Cdouble = Cdouble(0)
  r_drainage::Cdouble = Cdouble(0)
  r_root_decay::Cdouble = Cdouble(0)
  psi_min::Cdouble = Cdouble(0)
  alpha::Cdouble = Cdouble(0)
  f_soilwater::Cdouble = Cdouble(0)
  d_soil::NT10 = nzero(10)
  f_root::NT10 = nzero(10)
  dt::NT10 = nzero(10)
  thermal_cond::NT10 = nzero(10)
  theta_vfc::NT10 = nzero(10)
  theta_vwp::NT10 = nzero(10)
  fei::NT10 = nzero(10)
  Ksat::NT10 = nzero(10)
  psi_sat::NT10 = nzero(10)
  b::NT10 = nzero(10)
  density_soil::NT10 = nzero(10)
  f_org::NT10 = nzero(10)
  ice_ratio::NT10 = nzero(10)
  thetam::NT10 = nzero(10)
  thetam_prev::NT10 = nzero(10)
  temp_soil_p::NT10 = nzero(10)
  temp_soil_c::NT10 = nzero(10)
  f_ice::NT10 = nzero(10)
  psim::NT10 = nzero(10)
  thetab::NT10 = nzero(10)
  psib::NT10 = nzero(10)
  r_waterflow::NT10 = nzero(10)
  km::NT10 = nzero(10)
  Kb::NT10 = nzero(10)
  KK::NT10 = nzero(10)
  Cs::NT10 = nzero(10)
  lambda::NT10 = nzero(10)
  Ett::NT10 = nzero(10)
  G::NT10 = nzero(10)
end


@with_kw mutable struct ClimateData
  Srad::Cdouble = 0.0
  LR::Cdouble = 0.0
  temp::Cdouble = 0.0
  rh::Cdouble = 0.0
  rain::Cdouble = 0.0
  wind::Cdouble = 0.0
  dr_o::Cdouble = 0.0
  df_o::Cdouble = 0.0
  dr_u::Cdouble = 0.0
  df_u::Cdouble = 0.0
end


@with_kw mutable struct Results
  gpp_o_sunlit::Cdouble = 0.0
  gpp_u_sunlit::Cdouble = 0.0
  gpp_o_shaded::Cdouble = 0.0
  gpp_u_shaded::Cdouble = 0.0
  plant_resp::Cdouble = 0.0
  npp_o::Cdouble = 0.0
  npp_u::Cdouble = 0.0
  GPP::Cdouble = 0.0
  NPP::Cdouble = 0.0
  NEP::Cdouble = 0.0
  soil_resp::Cdouble = 0.0
  Net_Rad::Cdouble = 0.0

  SH::Cdouble = 0.0
  LH::Cdouble = 0.0
  Trans::Cdouble = 0.0
  Evap::Cdouble = 0.0
end


@with_kw mutable struct Cpools
  Ccd::NTuple{3,Cdouble} = nzero(3)
  Cssd::NTuple{3,Cdouble} = nzero(3)
  Csmd::NTuple{3,Cdouble} = nzero(3)
  Cfsd::NTuple{3,Cdouble} = nzero(3)
  Cfmd::NTuple{3,Cdouble} = nzero(3)
  Csm::NTuple{3,Cdouble} = nzero(3)
  Cm::NTuple{3,Cdouble} = nzero(3)
  Cs::NTuple{3,Cdouble} = nzero(3)
  Cp::NTuple{3,Cdouble} = nzero(3)
end


include("Constant.jl")
include("Leaf.jl")
include("OutputET.jl")

export Soil, ClimateData, Cpools, 
  Results, OutputET
export LeafRef
