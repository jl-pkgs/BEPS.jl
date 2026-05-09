export StateSeries, save_state!, split_vars
using OrderedCollections

const StateSeries{S<:NamedTuple,V<:NamedTuple} =
    NamedTuple{(:scalars, :vectors),Tuple{S,V}}


"""
预分配状态输出缓冲区
- scalar_fields: 标量字段名, e.g. (:z_water, :z_snow)
- vector_fields: 向量字段名, e.g. (:θ, :Tsoil_c)
"""
function StateSeries(::Val{SF}, ::Val{VF}, n_layer::Int, n_time::Int) where {SF,VF}
    Ns, Nv = length(SF), length(VF)
    s = NamedTuple{SF}(ntuple(_ -> Vector{Float64}(undef, n_time), Ns))
    v = NamedTuple{VF}(ntuple(_ -> Matrix{Float64}(undef, n_layer, n_time), Nv))
    (; scalars=s, vectors=v)
end

@generated function save_state!(out::StateSeries, st::StateBEPS, t::Int,
    ::Val{SF}, ::Val{VF}) where {SF,VF}
    nlayer = 5
    s_ex = [:(out.scalars.$(SF[i])[t] = st.$(SF[i])) for i in eachindex(SF)]
    v_ex = [:(copyto!(view(out.vectors.$(VF[i]), :, t), view(st.$(VF[i]), 1:$nlayer)))
            for i in eachindex(VF)]
    quote
        $(s_ex...)
        $(v_ex...)
        nothing
    end
end

save_state!(out::Nothing, st::StateBEPS, t::Int, ::Val{SF}, ::Val{VF}) where {SF,VF} = nothing


"""
    split_vars(vars) -> (sf::Val, vf::Val)

根据 VARS_SCALAR / VARS_VECTOR 将 vars 划分为标量和向量两组，并做合法性检查。
"""
function split_vars(vars)
    vars = Symbol.(vars)
    known = Set(fieldnames(StateBEPS))

    # 检查字段是否存在于 StateBEPS
    unknown = filter(v -> v ∉ known, vars)
    isempty(unknown) || error("字段不存在于 StateBEPS: $unknown")

    # 检查字段是否在可导出列表中
    exportable = Set((VARS_SCALAR..., VARS_VECTOR...))
    excluded = filter(v -> v ∉ exportable, vars)
    isempty(excluded) || @warn "字段在排除列表中，无法导出: $excluded"

    sf = Val(Tuple(v for v in vars if v ∈ Set(VARS_SCALAR)))
    vf = Val(Tuple(v for v in vars if v ∈ Set(VARS_VECTOR)))
    sf, vf
end

Base.Dict(nt::NamedTuple) = OrderedDict(pairs(nt))

function Base.getindex(out::StateSeries, i::Int)
    scalars = map(x -> x[i], out.scalars)
    vectors = map(x -> x[:, i], out.vectors)
    # Dict(:scalars => Dict(scalars), :vectors => Dict(vectors))
    (; scalars..., vectors...) |> Dict
end


function _print_section(io, nt, prefix)
    names = keys(nt)
    nw = maximum(length ∘ string, names)
    for (i, fname) in enumerate(names)
        c = i == length(names) ? "└" : "├"
        data = nt[fname]
        mb = round(Base.summarysize(data) / 1024^2, digits=2)
        print(io, "$(prefix)$(c)─ ")
        printstyled(io, rpad(string(fname), nw); bold=true, color=:blue)
        print(io, "  ", typeof(data), " | ", size(data), " | ")
        printstyled(io, "$mb Mb\n"; color=:blue, bold=true, underline=true)
    end
end

function Base.show(io::IO, ::MIME"text/plain", out::StateSeries)
    N = 52
    printstyled(io, "─"^N * "\n")
    # println(io, "─ scalars:")
    _print_section(io, out.scalars, "")
    printstyled(io, "─"^N * "\n")
    # println(io, "─ vectors:")
    _print_section(io, out.vectors, "")
    printstyled(io, "─"^N * "")
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", out::OrderedDict)
    # printstyled(io, "─"^N * "\n")
    # println(io, "─ scalars:")
    _print_section(io, out, "")
    # printstyled(io, "─"^N * "\n")
    # println(io, "─ vectors:")
    # _print_section(io, out.vectors, "")
    # printstyled(io, "─"^N * "")
    return nothing
end


export _print_section
