@with_kw mutable struct InterTempVars
  Tc_u::Vector{FT} = zeros(MAX_Loop)
  Ts0::Vector{FT} = zeros(MAX_Loop)
  Tsm0::Vector{FT} = zeros(MAX_Loop)
  Tsn0::Vector{FT} = zeros(MAX_Loop)
  Tsn1::Vector{FT} = zeros(MAX_Loop)
  Tsn2::Vector{FT} = zeros(MAX_Loop)
  Qhc_o::Vector{FT} = zeros(MAX_Loop)
  Qhc_u::Vector{FT} = zeros(MAX_Loop)
  Qhg::Vector{FT} = zeros(MAX_Loop)
  Wcl_o::Vector{FT} = zeros(MAX_Loop)
  Wcs_o::Vector{FT} = zeros(MAX_Loop) # the masses of liquid and snow on the canopy
  Xcl_o::Vector{FT} = zeros(MAX_Loop)
  Xcs_o::Vector{FT} = zeros(MAX_Loop) # the fraction of rain and snow on the canopy
  Wcl_u::Vector{FT} = zeros(MAX_Loop)
  Wcs_u::Vector{FT} = zeros(MAX_Loop) # the masses of rain and snow on the canopy
  Wg_snow::Vector{FT} = zeros(MAX_Loop)
  Xcl_u::Vector{FT} = zeros(MAX_Loop)
  Xcs_u::Vector{FT} = zeros(MAX_Loop) # the fraction of rain and snow on the canopy
  Ac_snow_o::Vector{FT} = zeros(MAX_Loop)
  Ac_snow_u::Vector{FT} = zeros(MAX_Loop)
  Xg_snow::Vector{FT} = zeros(MAX_Loop)
  rho_snow::Vector{FT} = zeros(MAX_Loop)
  alpha_v_sw::Vector{FT} = zeros(MAX_Loop)
  alpha_n_sw::Vector{FT} = zeros(MAX_Loop)
  r_rain_g::Vector{FT} = zeros(MAX_Loop)
  Trans_o::Vector{FT} = zeros(MAX_Loop)
  Trans_u::Vector{FT} = zeros(MAX_Loop)
  Eil_o::Vector{FT} = zeros(MAX_Loop)
  Eil_u::Vector{FT} = zeros(MAX_Loop)
  EiS_o::Vector{FT} = zeros(MAX_Loop)
  EiS_u::Vector{FT} = zeros(MAX_Loop)
  Evap_soil::Vector{FT} = zeros(MAX_Loop)
  Evap_SW::Vector{FT} = zeros(MAX_Loop)
  Evap_SS::Vector{FT} = zeros(MAX_Loop)
  lambda_snow::Vector{FT} = zeros(MAX_Loop)

  # 记录土壤温度
  
  Cs::Matrix{FT} = zeros(layer + 2, MAX_Loop)
  Tm::Matrix{FT} = zeros(layer + 2, MAX_Loop)
  G::Matrix{FT} = zeros(layer + 2, MAX_Loop)
end

function init_vars!(x::InterTempVars)
  # x.Tc_u .= 0.
  # x.Ts0 .= 0.
  # x.Tsn0 .= 0.
  # x.Tsm0 .= 0.
  # x.Tsn1 .= 0.
  # x.Tsn2 .= 0.
  # x.Qhc_o .= 0.
  # x.Qhc_u .= 0.
  # x.Qhg .= 0.
  x.Wcl_o .= 0.
  x.Wcs_o .= 0.
  x.Xcl_o .= 0.
  x.Xcs_o .= 0.
  x.Wcl_u .= 0.
  x.Wcs_u .= 0.
  x.Wg_snow .= 0.
  x.Xcl_u .= 0.
  x.Xcs_u .= 0.
  x.Xg_snow .= 0.0
  x.rho_snow .= 0.0

  # x.Ac_snow_o .= 0.
  # x.Ac_snow_u .= 0.  
  
  # x.alpha_v_sw .= 0.
  # x.alpha_n_sw .= 0.
  # x.r_rain_g .= 0.

  # x.Trans_o .= 0.
  # x.Trans_u .= 0.
  # x.Eil_o .= 0.
  # x.Eil_u .= 0.
  # x.EiS_o .= 0.
  # x.EiS_u .= 0.
  # x.Evap_soil .= 0.
  # x.Evap_SW .= 0.
  # x.Evap_SS .= 0.
  # x.lambda_snow .= 0.
  # x.Cs .= 0.
  # x.Tm .= 0.
  # x.G .= 0.
end
