import Parameters: @with_kw, @with_kw_noshow


const FW_VERSION = 1

const MAX_LAYERS = 10

const DEPTH_F = 6

const NOERROR = 0

const ERROR = 1

const PI = 3.1415926

const zero = 1.0e-10

const l_sta = 105

const l_end = 105

const p_sta = 101

const p_end = 101

const RTIMES = 24

const step = 3600.0

const kstep = 360.0

const kloop = 10

const MAX_Loop = 11

const layer = 5

const depth_f = 6

const CO2_air = 380.0    # ppm

const rho_a = 1.292    # density of air, kg m-3

const rho_w = 1025.0;  # density of water, kg m-3

const Cpd = 1004.65;

const Lv_solid = 2.83 * 1000000;  # the latent heat of evaporation from solid (snow/ice) at air temperature=Ta, in j+kkk/kg

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

dbl() = Cdouble(0)

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


@with_kw mutable struct Leaf
  o_sunlit::Cdouble = 0.0
  o_shaded::Cdouble = 0.0
  u_sunlit::Cdouble = 0.0
  u_shaded::Cdouble = 0.0
end

@with_kw mutable struct LeafRef{T}
  o_sunlit::T = Ref(0.0)
  o_shaded::T = Ref(0.0)
  u_sunlit::T = Ref(0.0)
  u_shaded::T = Ref(0.0)
end

function multiply!(Z::Leaf, X::Leaf, Y::Leaf)
  Z.o_sunlit = X.o_sunlit * Y.o_sunlit
  Z.o_shaded = X.o_shaded * Y.o_shaded
  Z.u_sunlit = X.u_sunlit * Y.u_sunlit
  Z.u_shaded = X.u_shaded * Y.u_shaded
end

function multiply!(Z::Leaf, X::LeafRef, Y::Leaf)
  Z.o_sunlit = X.o_sunlit[] * Y.o_sunlit
  Z.o_shaded = X.o_shaded[] * Y.o_shaded
  Z.u_sunlit = X.u_sunlit[] * Y.u_sunlit
  Z.u_shaded = X.u_shaded[] * Y.u_shaded
end


function init_leaf_struct(x::Leaf, replacement::Leaf)
  ccall((:init_leaf_struct, libbeps), Cvoid, (Ptr{Leaf}, Ptr{Leaf}), Ref(x), Ref(replacement))
end

function init_leaf_struct(x::Leaf, replacement::LeafRef)
  x.o_sunlit = replacement.o_sunlit[]
  x.o_shaded = replacement.o_shaded[]
  x.u_sunlit = replacement.u_sunlit[]
  x.u_shaded = replacement.u_shaded[]
end

function init_leaf_dbl(x::LeafRef, replacement::Float64)
  x.o_sunlit[] = replacement
  x.o_shaded[] = replacement
  x.u_sunlit[] = replacement
  x.u_shaded[] = replacement
end

function init_leaf_dbl(x::Leaf, replacement::Float64)
  x.o_sunlit = replacement
  x.o_shaded = replacement
  x.u_sunlit = replacement
  x.u_shaded = replacement
  # ccall((:init_leaf_dbl, libbeps), Cvoid, (Ptr{Leaf}, Cdouble), Ref(x), replacement)
end

function init_leaf_dbl2(x, overstory, understory)
  ccall((:init_leaf_dbl2, libbeps), Cvoid, (Ptr{Leaf}, Cdouble, Cdouble), Ref(x), overstory, understory)
end


include("DataType_ET.jl")

export Soil, ClimateData, Cpools, 
  Results, OutputET
export LeafRef
