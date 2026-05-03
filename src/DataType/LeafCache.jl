@with_kw mutable struct LeafCache
  init::Float64 = 0.0
  pc::PhotoConsts{Float64} = PhotoConsts(10.0 + 273.15) # 默认10°, 计算光合常量
  ac::AeroConsts{Float64} = AeroConsts()
  Ra::Radiation = Radiation()
  Cc_new::Leaf = Leaf(init)
  Cs_old::Leaf = Leaf(init)
  Cs_new::Leaf = Leaf(init)
  Ci_old::Leaf = Leaf(init)
  Tc_old::Leaf = Leaf(init)
  Tc_new::Leaf = Leaf(init)
  Gs_old::Leaf = Leaf(init)

  # to the reference height above the canopy
  Gc    ::Leaf = Leaf(init)  # total conductance for CO2 from the intercellular space of the leaves
  Gh    ::Leaf = Leaf(init)  # total conductance for heat transfer from the leaf surface
  Gw    ::Leaf = Leaf(init)  # total conductance for water from the intercellular space of the leaves
  Gww   ::Leaf = Leaf(init) # total conductance for water from the surface of the leaves

  Gs_new::Leaf = Leaf(init)
  Ac    ::Leaf = Leaf(init)
  Ci_new::Leaf = Leaf(init)

  Rn    ::Leaf = Leaf(init)
  Rns   ::Leaf = Leaf(init)
  Rnl   ::Leaf = Leaf(init)

  leleaf::Leaf = Leaf(init)
  GPP   ::Leaf = Leaf(init)
  LAI   ::Leaf = Leaf(init)
  PAI   ::Leaf = Leaf(init)
end

LeafCache(init) = LeafCache(; init)

# function reset!(l::LeafCache)
#   names = fieldnames(LeafCache)[2:end]
#   for name in names
#     x = getfield(l, name)
#     reset!(x)
#   end
# end
