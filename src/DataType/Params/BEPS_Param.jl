@bounds @with_kw_noshow mutable struct ParamBEPS{FT<:AbstractFloat}
  N::Int = 5
  dz::Vector{FT} = FT[0.05, 0.10, 0.20, 0.40, 1.25]  # 土壤层厚度 [m], BEPS V2023
  r_drainage::FT = Cdouble(0.50) | (0.2, 0.7)  # ? 地表排水速率（地表汇流），可考虑采用曼宁公式

  ψ_min::FT = Cdouble(33.0)  # [m], about 0.10~0.33 MPa开始胁迫点
  alpha::FT = Cdouble(0.4)   # [-], 土壤水限制因子参数，He 2017 JGR-B, Eq. 4

  hydraulic::ParamSoilHydraulicLayers{FT} = ParamSoilHydraulicLayers{FT,N}()
  thermal::ParamSoilThermalLayers{FT} = ParamSoilThermalLayers{FT,N}()

  veg::ParamVeg{FT} = ParamVeg{FT}()
end

# `kw...`: other params like, `r_drainage`
function ParamBEPS(VegType::Int, SoilType::Int; N::Int=5, FT=Float64, kw...)
  veg = InitParam_Veg(VegType; FT)
  hydraulic, thermal = InitParam_Soil(SoilType, N, FT)

  ψ_min = veg.is_bforest ? FT(10.0) : FT(33.0) # 开始胁迫点
  alpha = veg.is_bforest ? FT(1.5) : FT(0.4)   # 土壤水限制因子参数，He 2017 JGR-B, Eq. 4

  ParamBEPS{FT}(;
    N, kw..., ψ_min, alpha,
    hydraulic, thermal, veg
  )
end


# 这里应该加一个show function，打印模型参数信息
function Base.show(io::IO, model::M) where {M<:ParamBEPS}
  printstyled(io, "$M, N = $(model.N)\n", color=:blue, bold=true)

  fields_all = fieldnames(M)
  fields = setdiff(fields_all, [:N, :hydraulic, :thermal, :veg])

  n = length(fields)
  for i = 1:n
    field = fields[i]
    value = getfield(model, field)
    type = typeof(value)
    isa(value, Function) && (type = Function)
    println(io, "  $field\t: {$type} $value")
    # (i != n) && print(io, "\n")
  end

  ss = 60
  println(io, "-"^ss)
  printstyled(io, "Hydraulic: ", color=:blue, bold=true)
  print(io, model.hydraulic)

  println(io, "-"^ss)
  printstyled(io, "Thermal: ", color=:blue, bold=true)
  print(io, model.thermal)

  println(io, "-"^ss)
  printstyled(io, "Veg: ", color=:blue, bold=true)
  print(io, model.veg)
  print("-"^ss)
  return nothing
end
