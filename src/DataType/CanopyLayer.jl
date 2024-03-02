import Base.show

abstract type AbstractLayer{FT} end

"""
x = CanopyLayer{Float64}()
x = CanopyLayer{Ref{Float64}}()
"""
@with_kw_noshow mutable struct CanopyLayer{FT} <: AbstractLayer{FT}
  o::FT = FT(0.0) # overlayer
  u::FT = FT(0.0) # underlayer
  g::FT = FT(0.0) # ground
end

# 注意，如果是Ref，将共用相同的地址
CanopyLayer(o::FT) where {FT} = CanopyLayer{FT}(o, o, o)

CanopyLayer(o::FT, u::FT) where {FT} = CanopyLayer{FT}(; o, u)

CanopyLayer(x::CanopyLayer{FT}) = CanopyLayer{FT}(; x.o, x.u, x.g)

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


export CanopyLayer
