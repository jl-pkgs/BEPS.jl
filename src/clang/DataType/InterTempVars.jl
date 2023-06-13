@with_kw mutable struct InterTempVars
  Tc_u::Vector{FT} = zeros(MAX_Loop)
  Ts0::Vector{FT} = zeros(MAX_Loop)
  Tsn0::Vector{FT} = zeros(MAX_Loop)
  Tsm0::Vector{FT} = zeros(MAX_Loop)
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

  Cs::Matrix{FT} = zeros(layer + 2, MAX_Loop)
  Tm::Matrix{FT} = zeros(layer + 2, MAX_Loop)
  G::Matrix{FT} = zeros(layer + 2, MAX_Loop)
end
