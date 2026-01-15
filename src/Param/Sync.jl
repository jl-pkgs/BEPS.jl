# DBF or EBF, low constaint threshold
function Params2Soil!(soil::Soil, params::ParamBEPS{FT}; BF=false) where {FT}
  soil.ψ_min = BF ? 10.0 : 33.0 # [m], about 0.10~0.33 MPa开始胁迫点
  soil.alpha = BF ? 1.5 : 0.4   # He 2017 JGR-B, Eq. 4

  (; hydraulic, thermal, N) = params
  soil.n_layer = Cint(N)
  soil.dz[1:5] .= [0.05, 0.10, 0.20, 0.40, 1.25] # BEPS V2023, 土壤层厚度[m]

  soil.r_drainage = Cdouble(params.r_drainage)
  soil.r_root_decay = Cdouble(params.veg.r_root_decay)
  UpdateRootFraction!(soil) # 更新根系分布

  soil.ψ_min = Cdouble(params.ψ_min)
  soil.alpha = Cdouble(params.alpha)

  soil.θ_vfc[1:N] .= Cdouble.(hydraulic.θ_vfc)
  soil.θ_vwp[1:N] .= Cdouble.(hydraulic.θ_vwp)
  soil.θ_sat[1:N] .= Cdouble.(hydraulic.θ_sat)
  soil.K_sat[1:N] .= Cdouble.(hydraulic.K_sat)
  soil.ψ_sat[1:N] .= Cdouble.(hydraulic.ψ_sat)
  soil.b[1:N] .= Cdouble.(hydraulic.b)

  soil.κ_dry[1:N] .= Cdouble.(thermal.κ_dry)
  soil.ρ_soil[1:N] .= Cdouble.(thermal.ρ_soil)
  soil.V_SOM[1:N] .= Cdouble.(thermal.V_SOM)
end

Params2Soil!(soil::AbstractSoil, params::Nothing) = nothing

## 添加一个逆函数
function Soil2Params!(params::ParamBEPS{FT}, soil::Soil) where {FT}
  N = Int(soil.n_layer)
  params.N = N

  params.r_drainage = FT(soil.r_drainage)
  params.veg.r_root_decay = FT(soil.r_root_decay)
  params.ψ_min = FT(soil.ψ_min)
  params.alpha = FT(soil.alpha)

  if length(params.dz) != N; resize!(params.dz, N); end
  params.dz .= FT.(soil.dz[1:N])

  (; hydraulic, thermal) = params

  for field in (:θ_vfc, :θ_vwp, :θ_sat, :K_sat, :ψ_sat, :b)
    dest = getfield(hydraulic, field)
    src = getfield(soil, field)
    if length(dest) != N; resize!(dest, N); end
    dest .= FT.(src[1:N])
  end

  for field in (:κ_dry, :ρ_soil, :V_SOM)
    dest = getfield(thermal, field)
    src = getfield(soil, field)
    if length(dest) != N; resize!(dest, N); end
    dest .= FT.(src[1:N])
  end
  return params
end
