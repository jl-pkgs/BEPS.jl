export InterTempLeafs

@with_kw mutable struct InterTempLeafs
  x0::Float64 = 0.0
  Cc_new::Leaf = Leaf(x0)
  Cs_old::Leaf = Leaf(x0)
  Cs_new::Leaf = Leaf(x0)
  Ci_old::Leaf = Leaf(x0)
  Tc_old::Leaf = Leaf(x0)
  Tc_new::Leaf = Leaf(x0)
  Gs_old::Leaf = Leaf(x0)

  # to the reference height above the canopy
  Gc::Leaf = Leaf(x0)  # total conductance for CO2 from the intercellular space of the leaves
  Gh::Leaf = Leaf(x0)  # total conductance for heat transfer from the leaf surface 
  Gw::Leaf = Leaf(x0)  # total conductance for water from the intercellular space of the leaves
  Gww::Leaf = Leaf(x0) # total conductance for water from the surface of the leaves

  Gs_new::Leaf = Leaf(x0)
  Ac::Leaf = Leaf(x0)
  Ci_new::Leaf = Leaf(x0)

  Rn::Leaf = Leaf(x0)
  Rns::Leaf = Leaf(x0)
  Rnl::Leaf = Leaf(x0)

  leleaf::Leaf = Leaf(x0)
  GPP::Leaf = Leaf(x0)
  LAI::Leaf = Leaf(x0)
  PAI::Leaf = Leaf(x0)
end

InterTempLeafs(x0) = InterTempLeafs(; x0)

function reset!(l::InterTempLeafs)
  # reset!(l.Cc_new)
  # reset!(l.Cs_old)
  # reset!(l.Cs_new)
  # reset!(l.Ci_old)
  # reset!(l.Tc_old)
  # reset!(l.Tc_new)
  # reset!(l.Gs_old)
  # reset!(l.Gc)
  # reset!(l.Gh)
  # reset!(l.Gw)
  # reset!(l.Gww)
  # reset!(l.Gs_new)
  # reset!(l.Ac)
  # reset!(l.Ci_new)
  # reset!(l.Rn)
  # reset!(l.Rns)
  # reset!(l.Rnl)
  # reset!(l.leleaf)
  # reset!(l.GPP)
  # reset!(l.LAI)
  # reset!(l.PAI)

  # names = fieldnames(InterTempLeafs)[2:end]
  # for name in names
  #   x = getfield(l, name)
  #   reset(x)
  # end
end



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

  # Leafs
  TempLeafs::InterTempLeafs = InterTempLeafs(0.0)
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
