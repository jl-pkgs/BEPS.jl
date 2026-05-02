export StateSeries, save_state!
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
        printstyled(io, rpad(string(fname), nw); bold=true, color = :blue)
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
