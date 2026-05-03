# TODO: 可用StaticArrays改进
# they also abstract vector but with name
@with_kw mutable struct Leaf <: FieldVector{4,Cdouble}
  o_sunlit::Cdouble = 0.0
  o_shaded::Cdouble = 0.0
  u_sunlit::Cdouble = 0.0
  u_shaded::Cdouble = 0.0
end

Leaf(x::Cdouble) = Leaf(x,x,x,x)
Leaf(o::Cdouble, u::Cdouble) = Leaf(o, u)

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
