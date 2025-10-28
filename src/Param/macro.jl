export parameters, update!

using DataFrames

abstract type AbstractLayers{FT} end
abstract type AbstractModel{FT} end
abstract type AbstractBEPSmodel{FT} <: AbstractModel{FT} end

macro make_layers_struct(sname, sname_new=nothing)
  isnothing(sname_new) && (sname_new = Symbol(sname, :Layers))

  stype = getfield(__module__, sname)
  names_list = collect(fieldnames(stype))
  types_list = fieldtypes(stype)
  
  x = stype() # for default values
  values = map(fname -> getfield(x, fname), names_list)

  field_expressions = []
  # push!(field_expressions, :(ntime::Int = 100))
  for (i, fname) in enumerate(names_list)
    ftype = types_list[i]
    value = values[i]
    ftype <: AbstractVector && continue
    push!(field_expressions, :($fname::Vector{FT} = fill($value, N)))
  end

  quote
    @with_kw mutable struct $sname_new{FT,N} <: AbstractLayers{FT}
      $(field_expressions...)
    end
    has_definedbounds(x::$sname) = true
    has_definedbounds(x::$sname_new) = true

    function get_params(x::$sname_new{FT,N}; path=[]) where {FT,N}
      subtype = $sname
      # fields = fieldnames(Type)
      res = map(field -> begin
          value = getfield(x, field)
          _path = [path..., field]
          bound = bounds(subtype, field)
          map(i -> (; path=[_path..., i], name=field,
              value=value[i], type = eltype(value), bound=bound), 1:N)
        end, $names_list)
      vcat(res...)
    end
  end |> esc
end

has_definedbounds(x) = false


## 把 bounds 分解成字段路径和对应的约束
function split_bounds(x::S) where {S}
  function use_predef(field)
    # 如果是一个结构体，则采用递归的方式
    value = getfield(x, field)
    has_definedbounds(value) || isstructtype(typeof(value))
  end
  fields = fieldnames(S)
  (filter(use_predef, fields), filter(!use_predef, fields))
end


function get_params(x::S; path=[]) where {S}
  fs_predef, fs_macro = split_bounds(x)

  res_predef = map(field -> begin
      value = getfield(x, field)
      get_params(value; path=[path..., field])
    end, fs_predef)

  res_macro = map(field -> begin
      # @show bounds(x, field)
      value = getfield(x, field)
      (; path=[path..., field], name=field,
        value, type=eltype(value), bound=bounds(x, field))
    end, fs_macro)
  res = vcat(res_macro..., res_predef...)
  filter(x -> !isnothing(x.bound), res)
end


function update!(model::S, paths::Vector, values::Vector{FT},
  ; params::Union{Nothing,DataFrame}=nothing) where {S,FT}
  isnothing(params) && (params = Params(model))

  for (path, value) in zip(paths, values)
    rows = filter(row -> row.path == path, params)
    if isempty(rows)
      error("Parameter path $(path) not found in model!")
    elseif size(rows, 1) > 1
      error("Duplicated parameters found for path $(path)!")
    end
    update!(model, rows.path[1], value; type=rows.type[1])
  end
end

function update!(model::S, path::Vector, value::FT; type::Type) where {S,FT}
  if length(path) == 1
    # @show model, path[1], value
    setfield!(model, path[1], type(value))
  elseif length(path) > 1
    submodel = getfield(model, path[1]) # 

    if isa(submodel, Vector) # 如果是多模型
      models = submodel
      i = path[2]
      
      if typeof(models[i]) == FT
        models[i] = type(value)
        return
      end
      # 下面是应对Struct Vector
      update!(models[i], path[3:end], value; type)
    else
      update!(submodel, path[2:end], value; type)
    end
  end
end

parameters(model) = get_params(model) |> DataFrame
