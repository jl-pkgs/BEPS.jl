import DataFrames: DataFrame

function Base.sum(df::DataFrame)
  vals = [sum(df[!, c]) for c in names(df)]
  keys = names(df)
  list(keys, vals)
end

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
export readdlm
