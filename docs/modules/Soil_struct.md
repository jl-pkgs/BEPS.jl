# Soil 结构体使用报告

## 1. 定义位置

**`src/DataType/Soil.jl`**

```julia
@with_kw mutable struct Soil <: AbstractSoil
  # 状态变量: z_water, z_snow, θ, Tsoil_c, ice_ratio
  # 参数: dz, Ksat, ψ_sat, b, κ_dry, ρ_soil 等
end
```

继承自 `AbstractSoil`，使用 `Parameters.jl` 的 `@with_kw` 宏实现带默认值的初始化。

## 2. 使用 Soil 的模块清单

| 文件                                       | 作用                    |
| ------------------------------------------ | ----------------------- |
| `src/DataType/Soil.jl`                     | **定义** `Soil` 结构体  |
| `src/DataType/DataType.jl`                 | **导出** `Soil`         |
| `src/Param/Sync.jl`                        | `Sync_Param_to_Soil!`   |
| `src/SoilPhysics/SoilPhysics.jl`           | 土壤初始化和状态更新    |
| `src/SoilPhysics/UpdateHeatFlux.jl`        | 土壤热传导              |
| `src/SoilPhysics/UpdateSoilMoisture.jl`    | 土壤水分平衡            |
| `src/SoilPhysics/soil_water_factor_v2.jl`  | 水分胁迫因子计算        |
| `src/surface_temperature.jl`               | 地表能量平衡            |
| `src/inter_prg.jl`                         | 小时步积分循环          |
| `src/beps_main.jl`                         | **实例化** `Soil()`     |
| `src/beps_modern.jl`                       | 现代入口点              |
| `src/Param/deprecated/Init_Soil_Parameters.jl` | 初始化函数          |
| `test/modules/test-Soil.jl`                | 单元测试                |

## 3. 关键函数

| 函数                      | 文件                              | 作用               |
| ------------------------- | --------------------------------- | ------------------ |
| `Sync_Param_to_Soil!`     | `src/Param/Sync.jl`               | 同步模型参数到Soil |
| `Init_Soil_var`           | `deprecated/Init_Soil_Parameters.jl` | 综合初始化      |
| `UpdateSoilMoisture`      | `src/SoilPhysics/...`             | Darcy流和入渗      |
| `UpdateHeatFlux`          | `src/SoilPhysics/...`             | 层间热传导         |
| `surface_temperature!`    | `src/surface_temperature.jl`      | 地表-大气耦合      |
| `soil_water_factor_v2`    | `src/SoilPhysics/...`             | 植物水分胁迫       |
| `UpdateSoilThermalConductivity` | `src/SoilPhysics/...`       | 热导率更新         |
| `Update_Cs`               | `src/SoilPhysics/...`             | 体积热容更新       |

## 4. 生命周期

```
Soil() 实例化
    ↓
Init_Soil_var 初始化
    ↓
inter_prg 每小时更新循环:
    1. soil_water_factor_v2(soil)        → 计算水分胁迫
    2. UpdateSoilThermalConductivity(soil) → 更新热导率
    3. Update_Cs(soil)                   → 更新热容
    4. surface_temperature!(soil, ...)   → 地表能量平衡
    5. UpdateHeatFlux(soil, ...)         → 层间热传导
    6. UpdateSoilMoisture(soil, ...)     → 水分平衡
    ↓
输出结果
```

## 5. 字段分类

### 状态变量 [state]
- `z_water`, `z_snow` - 地表积水/积雪深度
- `θ[1:10]` - 土壤含水量
- `Tsoil_c[1:10]` - 当前土壤温度
- `ice_ratio[1:10]` - 冰比例
- `ψ[1:10]` - 土壤水势
- `G[1:10]` - 土壤热通量

### 参数 [param]
- `dz[1:10]` - 土壤层厚度
- `θ_sat`, `θ_vfc`, `θ_vwp` - 饱和/田间持水量/凋萎点
- `Ksat`, `ψ_sat`, `b` - 水力参数
- `κ_dry`, `ρ_soil`, `V_SOM` - 热力参数
- `ψ_min`, `alpha` - 水分胁迫参数

### 派生变量 [derived]
- `km`, `KK` - 水力导度
- `κ`, `Cs` - 热导率和热容
- `Ett[1:10]` - 每层蒸散量
- `f_soilwater` - 总水分胁迫因子
