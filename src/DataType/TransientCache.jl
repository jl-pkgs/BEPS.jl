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


@with_kw mutable struct TransientCache
  Tc_u::Vector{FT} = zeros(MAX_Loop)
  T_ground::Vector{FT} = zeros(MAX_Loop)      # 地表温度
  T_surf_mix::Vector{FT} = zeros(MAX_Loop)    # 混合表面温度
  T_surf_snow::Vector{FT} = zeros(MAX_Loop)   # 雪表面温度
  T_snow_L1::Vector{FT} = zeros(MAX_Loop)     # 雪层1温度
  T_snow_L2::Vector{FT} = zeros(MAX_Loop)     # 雪层2温度

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
  κ_snow::Vector{FT} = zeros(MAX_Loop)

  # 土壤温度和热通量
  Cs::Matrix{FT} = zeros(layer + 2, MAX_Loop)      # 土壤体积热容
  T_soil::Matrix{FT} = zeros(layer + 2, MAX_Loop)  # 土壤层温度
  G::Matrix{FT} = zeros(layer + 2, MAX_Loop)       # 土壤层热通量

  # Leafs
  TempLeafs::InterTempLeafs = InterTempLeafs(0.0)
end

function init_vars!(x::TransientCache)
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
end
