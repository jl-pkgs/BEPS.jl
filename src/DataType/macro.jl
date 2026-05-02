abstract type AbstractState{FT} end
abstract type AbstractFlux end

abstract type AbstractSeries{FT} end
abstract type AbstractFluxSeries{FT} <: AbstractSeries{FT} end
abstract type AbstractStateSeries{FT} <: AbstractSeries{FT} end


# ── 通用写入：series[i] = struct，按字段名静态展开 ─────────────────────────────
# 仿 setindex! 模式：编译期取两侧字段交集 → 内联赋值
@generated function Base.setindex!(series::SR, r::S, i::Int) where {FT<:AbstractFloat,
  SR<:AbstractSeries{FT}, S<:Union{AbstractFlux, AbstractState{FT}}}
  fs = intersect(fieldnames(SR), fieldnames(S))
  assigns = [:(@inbounds getfield(series, $(QuoteNode(f)))[i] = getfield(r, $(QuoteNode(f)))) for f in fs]
  return Expr(:block, assigns..., :(series))
end

# ── Series 定义宏 ─────────────────────────────────────────────────────────────
# 用法:
#   @DefFluxSeries FluxBEPS                         # 自动按命名约定 → FluxSeriesBEPS
#   @DefFluxSeries FluxSeries = Flux           # 显式指定 Series 名
#   @DefFluxSeries ETSeries = ETFlux extra1 extra2

macro DefFluxSeries(arg, extra_fields...)
  state_name, series_name = _parse_series_arg(arg)
  _def_series(__module__, state_name, extra_fields...; series_name, super_type=AbstractFluxSeries)
end

macro DefStateSeries(arg, extra_fields...)
  state_name, series_name = _parse_series_arg(arg)
  _def_series(__module__, state_name, extra_fields...; series_name, super_type=AbstractStateSeries)
end

# 解析: `Foo`              → (Foo, nothing)   走自动命名约定
#       `SeriesName = Foo` → (Foo, SeriesName)
function _parse_series_arg(arg)
  if isa(arg, Expr) && arg.head === :(=)
    return arg.args[2], arg.args[1]
  else
    return arg, nothing
  end
end

function _def_series(mod, StateStruct, extra_fields...; super_type, series_name=nothing)
  state_name = StateStruct
  output_name = something(series_name,
    Symbol(replace(string(state_name), r"^(Flux|State)" => s"\1Series")))

  state_type = getfield(mod, state_name)
  fieldnames_list = collect(fieldnames(state_type))
  field_types = fieldtypes(state_type)

  field_expressions = []

  # 1. ntime 字段（必须第一个）
  push!(field_expressions, :(ntime::Int = 100))

  # 2. 从源结构复制所有字段 → 转为 Vector{FT}（自动跳过 Vector 字段）
  for (i, fname) in enumerate(fieldnames_list)
    ftype = field_types[i]
    ftype <: AbstractVector && continue
    push!(field_expressions, :($fname::Vector{FT} = zeros(FT, ntime)))
  end

  # 3. 额外字段
  for extra in extra_fields
    if isa(extra, Symbol)
      push!(field_expressions, :($extra::Vector{FT} = zeros(FT, ntime)))
    elseif isa(extra, Expr) && extra.head === :(=)
      push!(field_expressions, extra)
    else
      error("额外字段格式错误：$(extra)")
    end
  end

  # State Series 提供 getindex 反向重建 State 结构
  scalar_fields = [fieldnames_list[i] for i in eachindex(fieldnames_list)
                   if !(field_types[i] <: AbstractVector)]

  getindex_expr = if super_type == AbstractStateSeries
    kw_args = [Expr(:kw, f, :(series.$f[i])) for f in scalar_fields]
    quote
      function Base.getindex(series::$output_name{FT}, i::Int) where {FT}
        $state_name{FT}($(kw_args...))
      end
    end
  else
    :()
  end

  return esc(quote
    @with_kw mutable struct $output_name{FT} <: $super_type{FT}
      $(field_expressions...)
    end
    $getindex_expr
  end)
end


# ── 单时刻结构定义宏（State / Flux）─────────────────────────────────────────
macro DefState(Struct, fields)
  field_syms = [f.value for f in fields.args if isa(f, QuoteNode)]
  field_expressions = [:($f::FT = 0.0) for f in field_syms]
  return esc(quote
    @with_kw mutable struct $Struct{FT<:AbstractFloat} <: AbstractState{FT}
      $(field_expressions...)
    end
  end)
end

macro DefFlux(Struct, fields)
  field_syms = [f.value for f in fields.args if isa(f, QuoteNode)]
  field_expressions = [:($f::FT = 0.0) for f in field_syms]
  return esc(quote
    @with_kw mutable struct $Struct{FT<:AbstractFloat} <: AbstractFlux
      $(field_expressions...)
    end
  end)
end

export AbstractFlux, AbstractState, AbstractFluxSeries, AbstractStateSeries, AbstractSeries
export @DefFluxSeries, @DefStateSeries, @DefState, @DefFlux


function Base.getindex(s::T, r::AbstractUnitRange{Int}) where {FT, T<:AbstractSeries{FT}}
  T(; ntime=length(r), (f => getfield(s, f)[r] for f in fieldnames(T) if f !== :ntime)...)
end

function Base.Matrix(res::AbstractSeries{T}) where {T<:Real}
  TYPE = typeof(res)
  names = fieldnames(TYPE)[2:end] |> collect
  data = map(i -> getfield(res, i), names)
  data = cat(data..., dims=2)
  data
end

function DataFrame(res::AbstractSeries{T}) where {T<:Real}
  TYPE = typeof(res)
  names = fieldnames(TYPE)[2:end] |> collect
  data = Matrix(res)
  DataFrame(data, names)
end
