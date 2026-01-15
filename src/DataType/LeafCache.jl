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
