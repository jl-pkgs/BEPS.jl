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

# function reset!(l::InterTempLeafs)
#   names = fieldnames(InterTempLeafs)[2:end]
#   for name in names
#     x = getfield(l, name)
#     reset!(x)
#   end
# end

export SurfaceMass

@with_kw mutable struct SurfaceMass{FT<:AbstractFloat}
  ρ_snow::FT = 0.0
  prcp_g::FT = 0.0

  depth_water::FT = 0.0
  depth_snow::FT = 0.0

  perc_snow_o::FT = 0.0 # Xcs_o
  perc_snow_u::FT = 0.0 # Xcs_u
  perc_snow_g::FT = 0.0 # Xcs_g

  perc_water_o::FT = 0.0 # Xcl_o
  perc_water_u::FT = 0.0 # Xcl_u

  area_snow_o::FT = 0.0 # Ac_snow_o
  area_snow_u::FT = 0.0 # Ac_snow_u

  pre_mass_water_o::FT = 0.0 # Wcl_o
  pre_mass_water_u::FT = 0.0 # Wcl_u

  mass_water_o::FT = 0.0 # Wcl_o
  mass_water_u::FT = 0.0 # Wcl_u
  mass_water_g::FT = 0.0 # Wcl_g

  mass_snow_o::FT = 0.0 # Wcs_o
  mass_snow_u::FT = 0.0 # Wcs_u
  mass_snow_g::FT = 0.0 # Wcs_g

  albedo_v_snow::FT = 0.0
  albedo_n_snow::FT = 0.0
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

  # the masses of liquid and snow on the canopy
  Wcl_o::Vector{FT} = zeros(MAX_Loop)
  Wcl_u::Vector{FT} = zeros(MAX_Loop)

  Wcs_o::Vector{FT} = zeros(MAX_Loop)
  Wcs_u::Vector{FT} = zeros(MAX_Loop)
  Wcs_g::Vector{FT} = zeros(MAX_Loop)

  # the fraction of rain and snow on the canopy
  Xcl_o::Vector{FT} = zeros(MAX_Loop)
  Xcl_u::Vector{FT} = zeros(MAX_Loop)

  Xcs_o::Vector{FT} = zeros(MAX_Loop)
  Xcs_u::Vector{FT} = zeros(MAX_Loop)
  Xcs_g::Vector{FT} = zeros(MAX_Loop)

  # area
  Ac_snow_o::Vector{FT} = zeros(MAX_Loop)
  Ac_snow_u::Vector{FT} = zeros(MAX_Loop)

  rho_snow::Vector{FT} = zeros(MAX_Loop)
  r_rain_g::Vector{FT} = zeros(MAX_Loop)

  alpha_v_sw::Vector{FT} = zeros(MAX_Loop)
  alpha_n_sw::Vector{FT} = zeros(MAX_Loop)

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
  x.Wcl_o .= 0.0
  x.Wcl_u .= 0.0

  x.Wcs_o .= 0.0
  x.Wcs_u .= 0.0
  x.Wcs_g .= 0.0

  x.Xcl_o .= 0.0
  x.Xcl_u .= 0.0

  x.Xcs_o .= 0.0
  x.Xcs_u .= 0.0
  x.Xcs_g .= 0.0

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
