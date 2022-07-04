# using BEPS
## MutableNamedTuples
# using MutableNamedTuples
# using Ipaper
import MutableNamedTuples: MutableNamedTuple

const MNT = MutableNamedTuple

# list(keys::Vector{Symbol}, values) = (; zip(keys, values)...)
"""
    list(keys::Vector{Symbol}, values)
    list(keys::Vector{<:AbstractString}, values)
    
# Examples
```julia
list([:dw, :betaw, :swmax, :a, :c, :kh, :uh]
```
"""
list(; kw...) = MNT(; kw...)


list(keys::Vector{Symbol}, values) = MNT(; zip(keys, values)...)

list(keys::Vector{<:AbstractString}, values) = list(Symbol.(keys), values)

list(keys::Vector{Symbol}) = list(Symbol.(keys), zeros(length(keys)))

list(keys::Tuple, values) = MNT(; zip(keys, values)...)
list(keys::Tuple) = list(keys, zeros(length(keys)))

list(keys::Vector{<:AbstractString}) = list(Symbol.(keys))

to_list = list;

Base.names(x::MNT) = keys(x) |> collect

# function Base.values(x::MNT)
#   @show "hello"
#   getindex.(values(getfield(x, :nt))) |> collect
#   # values(x) |> collect
# end


function add(x::MutableNamedTuple, y::MutableNamedTuple)
  list([keys(x)..., keys(y)...],
    [values(x)..., values(y)...],)
end

function Base.:(==)(x::MNT, y::MNT)
  if length(x) != length(y)
    return false
  end
  keys(x) == keys(y) && values(x) == values(y)
end

export list, add, MNT

# a = MNT{(:o_sunlit, :o_shaded, :u_sunlit, :u_shaded), NTuple{4, Base.RefValue{Int64}}}


## second version
import LabelledArrays
import LabelledArrays: LVector, LArray, symnames

const LA = LabelledArrays
LA.LVector(keys::Vector{Symbol}, values) = LVector(; zip(keys, values)...)
LA.LVector(keys::Tuple, values) = LVector(; zip(keys, values)...)

LA.LVector(keys::Vector{Symbol}) = LVector(keys, zeros(length(keys)))
LA.LVector(keys::Tuple) = LVector(keys, zeros(length(keys)))


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


export LVector
