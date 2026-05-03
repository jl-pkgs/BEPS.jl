export LeafCache, CacheSeries, save_cache!, split_cache_vars
export ALL_VARS_CACHE

@with_kw mutable struct LeafCache
  init::Float64 = 0.0
  pc::PhotoConsts{Float64} = PhotoConsts(10.0 + 273.15) # 默认10°, 计算光合常量
  ac::AeroConsts{Float64} = AeroConsts()
  Ra::Radiation = Radiation()
  # Cc_new::Leaf = Leaf(init)
  Cs_old::Leaf = Leaf(init)
  Cs_new::Leaf = Leaf(init)
  Ci_old::Leaf = Leaf(init)
  Ci_new::Leaf = Leaf(init)

  Tc_old::Leaf = Leaf(init)
  Tc_new::Leaf = Leaf(init)
  Gs_old::Leaf = Leaf(init)
  Gs_new::Leaf = Leaf(init)  # H2O 气孔导度 : [leaf intercellular]  -> [leaf surface]

  # to the reference height above the canopy
  Gc::Leaf = Leaf(init)      # CO2 总导度   : [leaf intercellular] -> [canopy reference height]
  Gh::Leaf = Leaf(init)      # heat 总导度  : [leaf surface]       -> [canopy reference height]
  Gw::Leaf = Leaf(init)      # H2O 总导度   : [leaf intercellular] -> [canopy reference height]
  Gw_wet::Leaf = Leaf(init)  # H2O 边界层导度: [wet leaf surface]   -> [canopy reference height]

  Ac::Leaf = Leaf(init)

  Rn::Leaf = Leaf(init)
  Rns::Leaf = Leaf(init)
  Rnl::Leaf = Leaf(init)

  leleaf::Leaf = Leaf(init)
  GPP::Leaf = Leaf(init)
  LAI::Leaf = Leaf(init)
  PAI::Leaf = Leaf(init)
end

LeafCache(init) = LeafCache(; init)

const ALL_VARS_CACHE = Tuple(
  f for (f, T) in zip(fieldnames(LeafCache), fieldtypes(LeafCache))
  if T == Leaf
)

const DEFAULT_CACHE_EXPORT = [:Tc_new, :Gs_new]

# function reset!(l::LeafCache)
#   names = fieldnames(LeafCache)[2:end]
#   for name in names
#     x = getfield(l, name)
#     reset!(x)
#   end
# end

function CacheSeries(::Val{VAR}, n_layer::Int, n_time::Int) where {VAR}
  bad = filter(v -> v ∉ ALL_VARS_CACHE, VAR)
  isempty(bad) || error("CacheSeries only supports Leaf fields: $bad")

  Ns = length(VAR)
  NamedTuple{VAR}(ntuple(_ -> Matrix{Float64}(undef, n_time, n_layer), Ns))
end

@generated function save_cache!(out::NamedTuple{VAR}, st::LeafCache, t::Int, ::Val{VAR}) where {VAR}
  bad = filter(v -> v ∉ ALL_VARS_CACHE, VAR)
  isempty(bad) || error("save_cache! only supports Leaf fields: $bad")

  v_ex = [:(copyto!(view(out.$(VAR[i]), t, :), st.$(VAR[i]))) for i in eachindex(VAR)]
  quote
    $(v_ex...)
    nothing
  end
end

function split_cache_vars(vars)
  vars = Symbol.(vars)
  known = Set(fieldnames(LeafCache))
  unknown = filter(v -> v ∉ known, vars)
  isempty(unknown) || error("字段不存在于 LeafCache: $unknown")

  exportable = Set(ALL_VARS_CACHE)
  excluded = filter(v -> v ∉ exportable, vars)
  isempty(excluded) || @warn "字段无法导出为 CacheSeries: $excluded"

  Val(Tuple(v for v in vars if v ∈ exportable))
end
