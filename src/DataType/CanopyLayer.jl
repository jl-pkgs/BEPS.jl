import Base.show

abstract type AbstractLayer{FT} end


## Layer3
"""
x = Layer3{Float64}()
x = Layer3{Ref{Float64}}()
"""
@with_kw_noshow mutable struct Layer3{FT} <: AbstractLayer{FT}
  o::FT = FT(0.0) # overlayer
  u::FT = FT(0.0) # underlayer
  g::FT = FT(0.0) # ground
end

# 注意，如果是Ref，将共用相同的地址
Layer3(o::FT) where {FT} = Layer3{FT}(o, o, o)

Layer3(o::FT, u::FT) where {FT} = Layer3{FT}(; o, u)

Layer3(x::Layer3{FT}) = Layer3{FT}(; x.o, x.u, x.g)

set!(x::Layer3{FT}, replacement::Layer3{FT}) where {FT} = begin
  x.o = replacement.o
  x.u = replacement.u
  x.g = replacement.g
end

## Layer2
@with_kw_noshow mutable struct Layer2{FT} <: AbstractLayer{FT}
  o::FT = FT(0.0) # overlayer
  u::FT = FT(0.0) # underlayer
end

# 注意，如果是Ref，将共用相同的地址
Layer2(o::FT) where {FT} = Layer2{FT}(o, o)

Layer2(x::Layer2{FT}) = Layer2{FT}(; x.o, x.u)

set!(x::Layer2{FT}, replacement::Layer2{FT}) where {FT} = begin
  x.o = replacement.o
  x.u = replacement.u
end

# Base_ops = ((:Base, :+), (:Base, :-), (:Base, :*), (:Base, :/))

# for (m, f) in Base_ops
#   @eval begin
#     $m.$f(a::Leaf, b::Leaf) = begin
#       Leaf(
#         $m.$f(a.o_sunlit, b.o_sunlit),
#         $m.$f(a.o_shaded, b.o_shaded),
#         $m.$f(a.u_sunlit, b.u_sunlit),
#         $m.$f(a.u_shaded, b.u_shaded)
#       )
#     end
#   end
# end


function Base.show(io::IO, x::AbstractLayer{FT}) where {FT}
  println(io, typeof(x))
  names = fieldnames(typeof(x))
  for i in 1:nfields(x)
    name = names[i]
    value = getfield(x, name)
    T = typeof(value)
    isa(value, Ref) && (value = value[])
    println(io, "$name: $T $(value)")
  end
end


export Layer3, Layer2
