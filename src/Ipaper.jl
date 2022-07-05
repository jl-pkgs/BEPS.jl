# using BEPS
## MutableNamedTuples
# using MutableNamedTuples
# using Ipaper
import MutableNamedTuples: MutableNamedTuple

const MNT = MutableNamedTuple

# mnt(keys::Vector{Symbol}, values) = (; zip(keys, values)...)
"""
    mnt(keys::Vector{Symbol}, values)
    mnt(keys::Vector{<:AbstractString}, values)
    
# Examples
```julia
mnt([:dw, :betaw, :swmax, :a, :c, :kh, :uh]
```
"""
mnt(; kw...) = MNT(; kw...)


mnt(keys::Vector{Symbol}, values) = MNT(; zip(keys, values)...)

mnt(keys::Vector{<:AbstractString}, values) = mnt(Symbol.(keys), values)

mnt(keys::Vector{Symbol}) = mnt(Symbol.(keys), zeros(length(keys)))

mnt(keys::Tuple, values) = MNT(; zip(keys, values)...)
mnt(keys::Tuple) = mnt(keys, zeros(length(keys)))

mnt(keys::Vector{<:AbstractString}) = mnt(Symbol.(keys))

to_list = mnt;

Base.names(x::MNT) = keys(x) |> collect

# function Base.values(x::MNT)
#   @show "hello"
#   getindex.(values(getfield(x, :nt))) |> collect
#   # values(x) |> collect
# end


function add(x::MutableNamedTuple, y::MutableNamedTuple)
  mnt([keys(x)..., keys(y)...],
    [values(x)..., values(y)...],)
end

function Base.:(==)(x::MNT, y::MNT)
  if length(x) != length(y)
    return false
  end
  keys(x) == keys(y) && values(x) == values(y)
end

export mnt, add, MNT

# a = MNT{(:o_sunlit, :o_shaded, :u_sunlit, :u_shaded), NTuple{4, Base.RefValue{Int64}}}


## second version
import LabelledArrays
import LabelledArrays: LVector, LArray, symnames

const LA = LabelledArrays
LA.LVector(keys::Vector{Symbol}, values) = LVector(; zip(keys, values)...)
LA.LVector(keys::Vector{<:AbstractString}, values) = LVector(; zip(Symbol.(keys), values)...)
LA.LVector(keys::Tuple, values) = LVector(; zip(keys, values)...)

LA.LVector(keys::Vector{Symbol}) = LVector(keys, zeros(length(keys)))
LA.LVector(keys::Tuple) = LVector(keys, zeros(length(keys)))

list = LVector;

Base.names(x::LArray) = symnames(typeof(x))


Base.:+(x::AbstractVector, y::AbstractVector) = x .+ y
Base.:-(x::AbstractVector, y::AbstractVector) = x .- y
Base.:*(x::AbstractVector, y::AbstractVector) = x .* y
Base.:/(x::AbstractVector, y::AbstractVector) = x ./ y

Base.:+(x::AbstractVector{T}, y::T) where {T<:Real} = x .+ y
Base.:-(x::AbstractVector{T}, y::T) where {T<:Real} = x .- y
Base.:*(x::AbstractVector{T}, y::T) where {T<:Real} = x .* y
Base.:/(x::AbstractVector{T}, y::T) where {T<:Real} = x ./ y

Base.:+(x::AbstractVector{T1}, y::T2) where {T1<:Real,T2<:Real} = x .+ y
Base.:-(x::AbstractVector{T1}, y::T2) where {T1<:Real,T2<:Real} = x .- y
Base.:*(x::AbstractVector{T1}, y::T2) where {T1<:Real,T2<:Real} = x .* y
Base.:/(x::AbstractVector{T1}, y::T2) where {T1<:Real,T2<:Real} = x ./ y

Base.:+(x::T, y::AbstractVector{T}) where {T<:Real} = x .+ y
Base.:-(x::T, y::AbstractVector{T}) where {T<:Real} = x .- y
Base.:*(x::T, y::AbstractVector{T}) where {T<:Real} = x .* y
Base.:/(x::T, y::AbstractVector{T}) where {T<:Real} = x ./ y

Base.:+(x::T1, y::AbstractVector{T2}) where {T1<:Real,T2<:Real} = x .+ y
Base.:-(x::T1, y::AbstractVector{T2}) where {T1<:Real,T2<:Real} = x .- y
Base.:*(x::T1, y::AbstractVector{T2}) where {T1<:Real,T2<:Real} = x .* y
Base.:/(x::T1, y::AbstractVector{T2}) where {T1<:Real,T2<:Real} = x ./ y


export LVector, list
