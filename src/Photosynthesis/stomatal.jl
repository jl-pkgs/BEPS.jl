using UnPack

"""
    ball_berry_gs(An::FT, RH::FT, cs::FT, params::ParamPhoto_Farquhar{FT}) where {FT<:AbstractFloat}

Ball-Berry 气孔导度模型

# Arguments
- `An`: 净光合速率 [μmol m-2 s-1]
- `RH`: 相对湿度 [0-1]
- `cs`: 叶表面 CO2 浓度 [μmol mol-1]
- `params`: 光合作用参数

# Returns
- `Gs`: 气孔导度 [mol m-2 s-1]

# References
- Ball et al., 1987
"""
@fastmath function ball_berry_gs(An::FT, RH::FT, cs::FT, params) where {FT<:AbstractFloat}
  @unpack g0_w, g1_w = params
  
  # Ball-Berry 模型: gs = g0 + g1 * An * hs / cs
  hs = RH  # 相对湿度
  gs = g0_w + g1_w * An * hs / cs
  
  return max(gs, g0_w)  # 确保不低于最小值
end

"""
    update_Gc(Gs::FT, gb_w::FT) where {FT<:AbstractFloat}

计算 CO2 总导度

# Arguments
- `Gs`: 气孔导度 [mol m-2 s-1]
- `gb_w`: 边界层导度 [mol m-2 s-1]

# Returns
- `Gc`: CO2 总导度 [mol m-2 s-1]

# Notes
- 1.4 是 CO2 与水汽的扩散系数比
- 1.37 是边界层扩散系数比
"""
@fastmath function update_Gc(Gs::FT, gb_w::FT) where {FT<:AbstractFloat}
  return 1.0 / (1.4 / Gs + 1.37 / gb_w)
end
