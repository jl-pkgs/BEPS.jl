import Base.show

## Layer3
"""
x = Layer3{Float64}()
x = Layer3{Ref{Float64}}()
"""
@with_kw_noshow mutable struct Layer3{FT} <: FieldVector{3,FT}
  o::FT = FT(0.0) # overlayer
  u::FT = FT(0.0) # underlayer
  g::FT = FT(0.0) # ground
end

## Layer2
@with_kw_noshow mutable struct Layer2{FT} <: FieldVector{2,FT}
  o::FT = FT(0.0) # overlayer
  u::FT = FT(0.0) # underlayer
end

const AbstractLayer{FT} = Union{Layer2{FT}, Layer3{FT}}


# 注意，如果是Ref，将共用相同的地址
Layer3(o::FT) where {FT} = Layer3{FT}(o, o, o)

Layer3(o::FT, u::FT) where {FT} = Layer3{FT}(; o, u, g = 0.0)

# 注意，如果是Ref，将共用相同的地址
Layer2(o::FT) where {FT} = Layer2{FT}(o, o)


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
