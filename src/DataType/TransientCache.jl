export LeafCache

@with_kw mutable struct LeafCache
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

LeafCache(x0) = LeafCache(; x0)

# function reset!(l::LeafCache)
#   names = fieldnames(LeafCache)[2:end]
#   for name in names
#     x = getfield(l, name)
#     reset!(x)
#   end
# end


@with_kw mutable struct TransientCache
  # 温度状态变量（需要历史访问 k-1）
  Tc_u::Vector{FT} = zeros(MAX_Loop)          # 下层冠层温度
  T_ground::Vector{FT} = zeros(MAX_Loop)      # 地表温度
  T_surf_mix::Vector{FT} = zeros(MAX_Loop)    # 混合表面温度
  T_surf_snow::Vector{FT} = zeros(MAX_Loop)   # 雪表面温度
  T_snow_L1::Vector{FT} = zeros(MAX_Loop)     # 雪层1温度
  T_snow_L2::Vector{FT} = zeros(MAX_Loop)     # 雪层2温度

  # 土壤温度和热通量（多层×时间步）
  Cs::Vector{FT} = zeros(layer + 2)      # 土壤体积热容
  G::Vector{FT} = zeros(layer + 2)       # 土壤层热通量

  # 叶片缓存（能量平衡迭代状态）
  leaf_cache::LeafCache = LeafCache(0.0)
end

function init_cache!(x::TransientCache)
  # 所有向量字段在创建时已经初始化为0，无需额外清零
  nothing
end
