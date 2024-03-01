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



function Base.sum(df::DataFrame)
  vals = [sum(df[!, c]) for c in names(df)]
  keys = names(df)
  list(keys, vals)
end

# for test
function nanmaximum(x::AbstractVector)
  inds = .!isnan.(x) |> findall
  !isempty(inds) ? maximum(x[inds]) : NaN
end

function Base.max(df::DataFrame)
  vals = [nanmaximum(df[!, c]) for c in names(df)]
  keys = names(df)
  list(keys, vals)
end

# import CSV
# fread(f) = DataFrame(CSV.File(f))
# fwrite(df, file) = begin
#   # dirname(file) |> check_dir
#   CSV.write(file, df)
# end
# export fread, fwrite, 
export LA
