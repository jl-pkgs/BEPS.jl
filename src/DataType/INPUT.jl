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

ClimateData(Srad, LR, temp, rh, rain, wind) = 
  ClimateData(; Srad, LR, temp, rh, rain, wind)


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
