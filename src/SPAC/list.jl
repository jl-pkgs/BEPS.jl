# export LVector, list
# export LA

# ## second version
# import LabelledArrays
# import LabelledArrays: LVector, LArray, symnames
# using DataFrames

# const LA = LabelledArrays
# LA.LVector(keys::Vector{Symbol}, values) = LVector(; zip(keys, values)...)
# LA.LVector(keys::Vector{<:AbstractString}, values) = LVector(; zip(Symbol.(keys), values)...)
# LA.LVector(keys::Tuple, values) = LVector(; zip(keys, values)...)

# LA.LVector(keys::Vector{Symbol}) = LVector(keys, zeros(length(keys)))
# LA.LVector(keys::Vector{<:AbstractString}) = LVector(Symbol.(keys), zeros(length(keys)))
# LA.LVector(keys::Tuple) = LVector(keys, zeros(length(keys)))

# list = LVector;
# Base.names(x::LArray) = symnames(typeof(x))


import DataFrames: AbstractDataFrame

function Base.sum(df::AbstractDataFrame)
  vals = [sum(df[!, c]) for c in names(df)]
  keys = names(df)
  NamedTuple{Tuple(Symbol.(keys))}(vals)
end


# for test
function nanmaximum(x::AbstractVector)
  inds = .!isnan.(x) |> findall
  !isempty(inds) ? maximum(x[inds]) : NaN
end

function Base.maximum(df::AbstractDataFrame)
  vals = [nanmaximum(df[!, c]) for c in names(df)]
  keys = names(df)
  NamedTuple{Tuple(Symbol.(keys))}(vals)
end

function _nanmaximum(x)
  x = collect(x)
  x = x[.!isnan.(x)]
  maximum(x)
end

export _nanmaximum

# import CSV
# fread(f) = DataFrame(CSV.File(f))
# fwrite(df, file) = begin
#   # dirname(file) |> check_dir
#   CSV.write(file, df)
# end
# export fread, fwrite, 
