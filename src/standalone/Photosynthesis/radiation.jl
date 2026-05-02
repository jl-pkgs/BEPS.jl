"""
    lai2!(LAI::FT, clumping::FT) where {FT<:AbstractFloat}

计算阳生和阴生叶片 LAI

# Arguments
- `LAI`: 总叶面积指数
- `clumping`: 叶片聚集指数

# Returns
- `LAI_sunlit`: 阳生叶片 LAI
- `LAI_shaded`: 阴生叶片 LAI

# References
- Chen et al., 2007
"""
@fastmath function lai2!(LAI::FT, clumping::FT) where {FT<:AbstractFloat}
  # 消光系数
  k = FT(0.5) / clumping
  
  # 阳生叶片 LAI
  LAI_sunlit = (1.0 - exp(-k * LAI)) / k
  
  # 阴生叶片 LAI
  LAI_shaded = LAI - LAI_sunlit
  
  return LAI_sunlit, LAI_shaded
end

"""
    calc_PAR_leaf(Srad::FT, LAI_layer::FT, k::FT=0.5) where {FT<:AbstractFloat}

计算叶片吸收的 PAR

# Arguments
- `Srad`: 入射短波辐射 [W m-2]
- `LAI_layer`: 累积 LAI 到该层
- `k`: 消光系数（默认 0.5）

# Returns
- `PAR`: 光合有效辐射 [μmol m-2 s-1]

# Notes
- 假设 PAR 占短波辐射的 50%
- 转换系数: 1 W m-2 ≈ 4.6 μmol m-2 s-1
"""
@fastmath function calc_PAR_leaf(Srad::FT, LAI_layer::FT, k::FT=FT(0.5)) where {FT<:AbstractFloat}
  # PAR 占短波辐射的 50%
  PAR_total = FT(0.5) * Srad * FT(4.6)  # [μmol m-2 s-1]
  
  # 该层吸收的 PAR（简化）
  PAR_leaf = PAR_total * exp(-k * LAI_layer)
  
  return PAR_leaf
end
