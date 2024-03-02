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


export @reset, reset!
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


function update_Gw!(Gw::Leaf, Gs_new::Leaf, Ga_o, Ga_u, Gb_o, Gb_u)
  Gw.o_sunlit = 1.0 / (1.0 / Ga_o + 1.0 / Gb_o + 1.0 / Gs_new.o_sunlit)  #conductance for water
  Gw.o_shaded = 1.0 / (1.0 / Ga_o + 1.0 / Gb_o + 1.0 / Gs_new.o_shaded)
  Gw.u_sunlit = 1.0 / (1.0 / Ga_u + 1.0 / Gb_u + 1.0 / Gs_new.u_sunlit)
  Gw.u_shaded = 1.0 / (1.0 / Ga_u + 1.0 / Gb_u + 1.0 / Gs_new.u_shaded)
end

function update_Gc!(Gc::Leaf, Gs_new::Leaf, Ga_o, Ga_u, Gb_o, Gb_u)
  Gc.o_sunlit = 1.0 / (1.0 / Ga_o + 1.4 / Gb_o + 1.6 / Gs_new.o_sunlit) # conductance for CO2
  Gc.o_shaded = 1.0 / (1.0 / Ga_o + 1.4 / Gb_o + 1.6 / Gs_new.o_shaded)
  Gc.u_sunlit = 1.0 / (1.0 / Ga_u + 1.4 / Gb_u + 1.6 / Gs_new.u_sunlit)
  Gc.u_shaded = 1.0 / (1.0 / Ga_u + 1.4 / Gb_u + 1.6 / Gs_new.u_shaded)
end
