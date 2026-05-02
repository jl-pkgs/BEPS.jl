using Parameters

"""
    LeafPhoto{FT}

简化的叶片光合作用变量（阳生/阴生）

# Fields
- `sunlit::FT`: 阳生叶片值
- `shaded::FT`: 阴生叶片值
"""
@with_kw mutable struct LeafPhoto{FT<:AbstractFloat}
  sunlit::FT = FT(0.0)
  shaded::FT = FT(0.0)
end

"""
    PhotoResult{FT}

光合作用计算结果

# Fields
- `Gc::LeafPhoto{FT}`: CO2 导度 [mol m-2 s-1]
- `An::LeafPhoto{FT}`: 净光合速率 [μmol m-2 s-1]
- `Ci::LeafPhoto{FT}`: 胞间 CO2 浓度 [μmol mol-1]
- `Gs::LeafPhoto{FT}`: 气孔导度 [mol m-2 s-1]
- `Rd::LeafPhoto{FT}`: 暗呼吸速率 [μmol m-2 s-1]
"""
@with_kw mutable struct PhotoResult{FT<:AbstractFloat}
  Gc::LeafPhoto{FT} = LeafPhoto{FT}()
  An::LeafPhoto{FT} = LeafPhoto{FT}()
  Ci::LeafPhoto{FT} = LeafPhoto{FT}()
  Gs::LeafPhoto{FT} = LeafPhoto{FT}()
  Rd::LeafPhoto{FT} = LeafPhoto{FT}()
end
