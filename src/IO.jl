
import DataFrames: DataFrame
import CSV
using DelimitedFiles: readdlm


function Base.sum(df::DataFrame)
  vals = [sum(df[!, c]) for c in names(df)]
  keys = names(df)
  list(keys, vals)
end


# import NaNStatistics: nanmaximum

function Base.max(df::DataFrame)
  vals = [maximum(df[!, c]) for c in names(df)]
  keys = names(df)
  list(keys, vals)
end

fread(f) = DataFrame(CSV.File(f))
fwrite(df, file) = begin
  # dirname(file) |> check_dir
  CSV.write(file, df)
end


export fread, fwrite, readdlm
