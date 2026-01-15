export LeafCache
export @reset, reset!

# TODO: 可用StaticArrays改进
# they also abstract vector but with name
@with_kw mutable struct Leaf
  o_sunlit::Cdouble = 0.0
  o_shaded::Cdouble = 0.0
  u_sunlit::Cdouble = 0.0
  u_shaded::Cdouble = 0.0
end

Leaf(x::Cdouble) = Leaf(x,x,x,x)
Leaf(o::Cdouble, u::Cdouble) = Leaf(o, u)


Base_ops = ((:Base, :+), (:Base, :-), (:Base, :*), (:Base, :/))

for (m, f) in Base_ops
  @eval begin
    $m.$f(a::Leaf, b::Leaf) = begin
      Leaf(
        $m.$f(a.o_sunlit, b.o_sunlit),
        $m.$f(a.o_shaded, b.o_shaded),
        $m.$f(a.u_sunlit, b.u_sunlit),
        $m.$f(a.u_shaded, b.u_shaded)
      )
    end
  end
end


macro reset(x)
  quote
    reset!($(esc(x)))
  end
end
reset!(x::Leaf) = init_leaf_dbl(x, 0.0)

function set!(x::Leaf, replacement::Leaf)
  x.o_sunlit = replacement.o_sunlit
  x.o_shaded = replacement.o_shaded
  x.u_sunlit = replacement.u_sunlit
  x.u_shaded = replacement.u_shaded
end

function init_leaf_dbl(x::Leaf, replacement::Float64)
  x.o_sunlit = replacement
  x.o_shaded = replacement
  x.u_sunlit = replacement
  x.u_shaded = replacement
end

function init_leaf_dbl2(x::Leaf, overstory, understory)
  x.o_sunlit = overstory
  x.o_shaded = overstory
  x.u_sunlit = understory
  x.u_shaded = understory
end


function multiply!(Z::Leaf, X::Leaf, Y::Leaf)
  Z.o_sunlit = X.o_sunlit * Y.o_sunlit
  Z.o_shaded = X.o_shaded * Y.o_shaded
  Z.u_sunlit = X.u_sunlit * Y.u_sunlit
  Z.u_shaded = X.u_shaded * Y.u_shaded
end



@with_kw mutable struct LeafCache
  init::Float64 = 0.0
  Cc_new::Leaf = Leaf(init)
  Cs_old::Leaf = Leaf(init)
  Cs_new::Leaf = Leaf(init)
  Ci_old::Leaf = Leaf(init)
  Tc_old::Leaf = Leaf(init)
  Tc_new::Leaf = Leaf(init)
  Gs_old::Leaf = Leaf(init)

  # to the reference height above the canopy
  Gc::Leaf = Leaf(init)  # total conductance for CO2 from the intercellular space of the leaves
  Gh::Leaf = Leaf(init)  # total conductance for heat transfer from the leaf surface 
  Gw::Leaf = Leaf(init)  # total conductance for water from the intercellular space of the leaves
  Gww::Leaf = Leaf(init) # total conductance for water from the surface of the leaves

  Gs_new::Leaf = Leaf(init)
  Ac::Leaf = Leaf(init)
  Ci_new::Leaf = Leaf(init)

  Rn::Leaf = Leaf(init)
  Rns::Leaf = Leaf(init)
  Rnl::Leaf = Leaf(init)

  leleaf::Leaf = Leaf(init)
  GPP::Leaf = Leaf(init)
  LAI::Leaf = Leaf(init)
  PAI::Leaf = Leaf(init)
end

LeafCache(init) = LeafCache(; init)

# function reset!(l::LeafCache)
#   names = fieldnames(LeafCache)[2:end]
#   for name in names
#     x = getfield(l, name)
#     reset!(x)
#   end
# end
