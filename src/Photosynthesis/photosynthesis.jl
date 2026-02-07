module Photosynthesis

using UnPack
using Parameters
# ParamPhoto_Farquhar 将在 BEPS 模块中定义

# 导入类型
include("types.jl")

# 导入子模块
include("helper.jl")
include("temperature.jl")
include("radiation.jl")
include("stomatal.jl")
include("core.jl")

# 导出公共 API
export LeafPhoto, PhotoResult
export photosynthesis

"""
    photosynthesis(Tair::FT, RH::FT, Srad::FT, LAI::FT, 
                   params; 
                   ca::FT=380.0, 
                   β_soil::FT=1.0,
                   gb_w::FT=0.01) where {FT<:AbstractFloat}

简化的光合作用模型（假设 T_leaf = T_air）

# Arguments
- `Tair`: 气温 [°C]
- `RH`: 相对湿度 [%]
- `Srad`: 短波辐射 [W m-2]
- `LAI`: 叶面积指数
- `params`: 光合作用参数（ParamPhoto_Farquhar）
- `ca`: 大气 CO2 浓度 [μmol mol-1] (默认 380)
- `β_soil`: 土壤水分胁迫因子 [0-1] (默认 1.0，无胁迫)
- `gb_w`: 边界层导度 [mol m-2 s-1] (默认 0.01)

# Returns
- `PhotoResult`: 包含 An, Gs, Ci, Gc, Rd

# Examples
```julia
using BEPS

params = InitParam_Photo_Farquhar()
result = photosynthesis(25.0, 60.0, 500.0, 2.0, params)
println("阳生叶片 An: ", result.An.sunlit, " μmol m-2 s-1")
```

# References
- Farquhar et al., 1980
- Ball et al., 1987
- Chen et al., 2007
"""
function photosynthesis(Tair::FT, RH::FT, Srad::FT, LAI::FT,
                       params;
                       ca::FT=FT(380.0),
                       β_soil::FT=FT(1.0),
                       gb_w::FT=FT(0.01)) where {FT<:AbstractFloat}
  # 计算气象变量
  ea = cal_ea(Tair, RH)
  rho_a = cal_rho_a(Tair, ea)
  T = Tair + FT(273.15)
  
  # LAI 分配到阳生/阴生叶片
  LAI_sunlit, LAI_shaded = lai2!(LAI, params.clumping)
  
  # 计算叶片吸收的 PAR
  PAR_sunlit = calc_PAR_leaf(Srad, FT(0.0))
  PAR_shaded = calc_PAR_leaf(Srad, LAI_sunlit)
  
  result = PhotoResult{FT}()
  
  # 迭代求解阳生叶片
  MAX_ITERATIONS = 15
  CONVERGENCE_THRESHOLD = FT(1.0)
  
  ci = ca * FT(0.7)
  An = FT(0.0)
  Rd = FT(0.0)
  gs = FT(0.0)
  gc = FT(0.0)
  
  for iter in 1:MAX_ITERATIONS
    An, Rd = farquhar_model(T, PAR_sunlit, ci, params)
    An = An * β_soil
    cs = ca
    gs = ball_berry_gs(max(An, 0.0), RH / FT(100.0), cs, params)
    gc = update_Gc(gs, gb_w)
    ci_new = ca - An / gc
    
    if abs(ci_new - ci) < CONVERGENCE_THRESHOLD
      ci = ci_new
      break
    end
    ci = ci_new
  end
  
  result.An.sunlit = An
  result.Gs.sunlit = gs
  result.Ci.sunlit = ci
  result.Rd.sunlit = Rd
  result.Gc.sunlit = gc
  
  # 迭代求解阴生叶片
  ci = ca * FT(0.7)
  An = FT(0.0)
  Rd = FT(0.0)
  gs = FT(0.0)
  gc = FT(0.0)
  
  for iter in 1:MAX_ITERATIONS
    An, Rd = farquhar_model(T, PAR_shaded, ci, params)
    An = An * β_soil
    cs = ca
    gs = ball_berry_gs(max(An, 0.0), RH / FT(100.0), cs, params)
    gc = update_Gc(gs, gb_w)
    ci_new = ca - An / gc
    
    if abs(ci_new - ci) < CONVERGENCE_THRESHOLD
      ci = ci_new
      break
    end
    ci = ci_new
  end
  
  result.An.shaded = An
  result.Gs.shaded = gs
  result.Ci.shaded = ci
  result.Rd.shaded = Rd
  result.Gc.shaded = gc
  
  return result
end

end # module
