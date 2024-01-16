export LVector, list

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
