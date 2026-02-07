# BEPS.jl 参数管理 - AI 友好版

## 核心概念

BEPS.jl 参数系统有三个关键要素：

1. **默认值**: 参数的初始值
2. **边界**: 参数的允许范围 (最小值, 最大值)
3. **传递**: 默认值从定义到使用的流程

---

## 1. 参数定义语法

```julia
@bounds @with_kw mutable struct ParamVeg{FT<:AbstractFloat}
  LAI_max_o::FT = 4.5 | (0.1, 7.0)
end
```

**语法解析**:
- `@bounds`: 添加边界元数据的宏
- `@with_kw`: Parameters.jl 宏，支持关键字参数
- `FT`: 类型参数 (Float32/Float64)
- `4.5`: 默认值
- `| (0.1, 7.0)`: 边界约束

**等价于**:
```julia
LAI_max_o::FT = 4.5
bounds(ParamVeg, :LAI_max_o) = (0.1, 7.0)
```

---

## 2. 默认值的三种来源

### 来源 1: 硬编码 (Params.jl)
```julia
@bounds @with_kw mutable struct ParamVeg{FT<:AbstractFloat}
  LAI_max_o::FT = 4.5 | (0.1, 7.0)  # 硬编码默认值
end
```
- **用途**: 通用常数
- **优先级**: 低 (可被覆盖)

### 来源 2: JSON 文件 (ParamVeg.json)
```json
{
  "ENF": {"VCmax25": 62.5},
  "DBF": {"VCmax25": 57.7}
}
```
```julia
veg = InitParam_Veg(6)  # 加载 DBF 类型参数
# veg.VCmax25 = 57.7 (来自 JSON，覆盖硬编码)
```
- **用途**: 植被类型特定参数
- **优先级**: 中

### 来源 3: 数据表 (GlobalData.jl)
```julia
const SOIL_PARAMS = [
  (name="sand", θ_sat=0.437, b=[1.7, 1.9, ...]),
  (name="loam", θ_sat=0.463, b=[4.5, 4.7, ...])
]

hydraulic, thermal = InitParam_Soil(4, 5, Float64)  # Loam, 5层
```
- **用途**: 土壤类型特定参数
- **优先级**: 中

---

## 3. 参数传递流程

```
┌─────────────────────────────────────────────────────────────┐
│ Step 1: 定义默认值                                          │
├─────────────────────────────────────────────────────────────┤
│ 硬编码: Params.jl          → LAI_max_o = 4.5               │
│ JSON: ParamVeg.json       → VCmax25 = 57.7 (DBF)           │
│ 数据表: GlobalData.jl     → θ_sat = 0.463 (Loam)           │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ Step 2: 初始化参数                                          │
├─────────────────────────────────────────────────────────────┤
│ InitParam_Veg(6)         → ParamVeg{FT}                    │
│ InitParam_Soil(4, 5)     → ParamSoilHydraulicLayers{FT,5}  │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ Step 3: 创建模型                                            │
├─────────────────────────────────────────────────────────────┤
│ model = ParamBEPS(VegType=6, SoilType=4; N=5)              │
│   ├─ veg::ParamVeg{FT}                                     │
│   ├─ hydraulic::ParamSoilHydraulicLayers{FT,5}             │
│   └─ thermal::ParamSoilThermalLayers{FT,5}                 │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ Step 4: 提取参数                                            │
├─────────────────────────────────────────────────────────────┤
│ df = parameters(model)                                      │
│ → DataFrame{path, name, value, type, bound}                │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ Step 5: 更新参数                                            │
├─────────────────────────────────────────────────────────────┤
│ update!(model, [:veg, :VCmax25], 100.0)                    │
│ → model.veg.VCmax25 = 100.0                                │
└─────────────────────────────────────────────────────────────┘
```

---

## 4. 常用操作

### 创建模型
```julia
model = ParamBEPS(VegType=6, SoilType=4; N=5)
# VegType: 6=DBF (落叶阔叶林)
# SoilType: 4=Loam (壤土)
# N: 5 层土壤
```

### 查看参数
```julia
# 方法 1: 打印模型
display(model)

# 方法 2: 提取为 DataFrame
df = parameters(model)
# df 包含: path, name, value, type, bound
```

### 更新参数
```julia
# 单个参数
update!(model, [:veg, :VCmax25], 100.0)

# 批量更新
paths = [[:veg, :VCmax25], [:hydraulic, :b, 3]]
values = [100.0, 4.5]
update!(model, paths, values; params=df)
```

### 参数优化
```julia
model_opt, rmse = beps_optimize(d, lai, model, obs_ET;
  col_sim=:ET, maxn=2000)
```

---

## 5. 植被和土壤类型

### 植被类型 (VegType)
```julia
1  = ENF  # 常绿针叶林
6  = DBF  # 落叶阔叶林
9  = EBF  # 常绿阔叶林
13 = Shrub # 灌木
40 = C4    # C4 草地
```

### 土壤类型 (SoilType)
```julia
1  = sand         # 沙土
4  = loam         # 壤土
11 = clay         # 黏土
```

---

## 6. 添加新参数的步骤

假设要添加新参数 `leaf_thickness`:

### Step 1: 在结构体中定义
```julia
# 文件: src/DataType/Params/Params.jl
@bounds @with_kw mutable struct ParamVeg{FT<:AbstractFloat}
  leaf_thickness::FT = 0.2 | (0.1, 0.5)  # [mm]
end
```

### Step 2: 在 JSON 中添加值
```json
// 文件: src/DataType/Params/data/ParamVeg.json
{
  "ENF": {"leaf_thickness": 0.25},
  "DBF": {"leaf_thickness": 0.18}
}
```

### Step 3: 在 InitParam_Veg 中加载
```julia
# 文件: src/DataType/Params/Param_Init.jl
function InitParam_Veg(lc::Int=1; FT=Float64)
  v = veg_data[type_str]
  return ParamVeg{FT}(
    leaf_thickness = FT(v["leaf_thickness"])
  )
end
```

### Step 4: 测试
```julia
veg = InitParam_Veg(1)
@test veg.leaf_thickness == 0.25
```

---

## 7. 关键要点总结

1. **参数定义**: `字段名::类型 = 默认值 | (最小值, 最大值)`
2. **默认值来源**: 硬编码、JSON、数据表
3. **传递流程**: 定义 → 初始化 → 模型 → 提取 → 更新
4. **边界查询**: `bounds(ParamVeg{Float64}, :LAI_max_o)`
5. **参数提取**: `df = parameters(model)`
6. **参数更新**: `update!(model, path, value)`

---

**这就是 BEPS.jl 参数管理的全部核心内容。**
