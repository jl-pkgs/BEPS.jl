# BEPS.jl 代码风格指南

> 基于 Linux 极简主义和 Julia 最佳实践

## 核心原则

1. **极简主义**: 代码简洁、清晰、无冗余
2. **类型安全**: 充分利用 Julia 类型系统
3. **性能优先**: 使用 `@inbounds`、`@fastmath` 等优化
4. **文档完善**: 英文文档字符串，中文行内注释
5. **物理严谨**: 数值截断、单位标注、边界检查

---

## 1. 命名规范

### 1.1 变量命名

```julia
# 希腊字母用于物理量
θ         # 土壤含水量 [m³/m³]
ψ         # 土壤基质势 [m]
κ         # 热导率 [W/m/K]
ρ         # 密度 [kg/m³]

# 后缀约定
Tsoil_p   # 上一时间步 (previous)
Tsoil_c   # 当前时间步 (current)
θ_sat     # 饱和值 (saturation)
θ_vwp     # 萎蔫点 (volumetric wilting point)
θ_vfc     # 田间持水量 (volumetric field capacity)

# 前缀约定
f_soilwater  # 分数/因子 (fraction)
r_drainage   # 速率 (rate)
z_snow       # 深度/厚度 (depth)

# 缩写词
LAI      # Leaf Area Index
VCmax    # Maximum carboxylation rate
ET       # Evapotranspiration
GPP      # Gross Primary Productivity
```

### 1.2 函数命名

```julia
# 动词开头，描述功能
InitParam_Veg()           # 初始化植被参数
UpdateSoilMoisture()      # 更新土壤水分
cal_rho_a()               # 计算空气密度 (cal = calculate)
solve_canopy_energy_balance!()  # 求解冠层能量平衡

# 后缀约定
*_jl                      # Julia 实现（区别于 C 版本）
!                        # 原地修改函数
```

### 1.3 类型命名

```julia
# PascalCase，描述性名称
ParamVeg                  # 植被参数
ParamSoilHydraulic        # 土壤水力参数
StateBEPS                 # BEPS 状态
Met                       # 气象输入

# 类型参数
FT<:AbstractFloat         # 浮点类型参数
T<:Cdouble                # 双精度浮点
S<:Union{StateBEPS,Soil}  # 联合类型
```

---

## 2. 结构体定义

### 2.1 使用 @with_kw 宏

```julia
using Parameters
import FieldMetadata: @metadata

@metadata bounds nothing

@bounds @with_kw mutable struct ParamVeg{FT<:AbstractFloat}
  # 字段::类型 = 默认值 | (最小值, 最大值)
  LAI_max_o::FT = 4.5 | (0.1, 7.0)
  VCmax25::FT = 89.45 | (5.0, 200.0)

  # 布尔值
  has_understory::Bool = true
  is_bforest::Bool = false
end
```

### 2.2 字段组织顺序

```julia
@bounds @with_kw mutable struct StateBEPS
  # 1. 维度和层数
  n_layer::Cint = 5
  dz::Vector{Float64} = zeros(10)

  # 2. 温度状态
  Tsnow_c::SnowLand{FT}
  Tsoil_c::Vector{Float64} = zeros(10)
  Tsoil_p::Vector{Float64} = zeros(10)

  # 3. 水分状态
  θ::Vector{Float64} = zeros(10)
  z_snow::Cdouble = 0
  z_water::Cdouble = 0

  # 4. 物理属性
  κ::Vector{Float64} = zeros(10)
  Cv::Vector{Float64} = zeros(10)

  # 5. 通量
  G::Vector{Float64} = zeros(10)
  Ett::Vector{Float64} = zeros(10)

  # 6. 胁迫因子
  f_soilwater::Cdouble = 0
  f_root::Vector{Float64} = zeros(10)
end
```

### 2.3 注释规范

```julia
@bounds @with_kw mutable struct ParamVeg{FT<:AbstractFloat}
  # 简短注释放在字段后面
  LAI_max_o::FT = 4.5 | (0.1, 7.0)  # LAI max for overstory

  # 复杂注释放在字段上方
  # decay_rate_of_root_distribution
  # 根系分布衰减率，控制根系在土壤中的垂直分布
  r_root_decay::FT = 0.95 | (0.85, 0.999)

  # 单位标注
  z_canopy_o::FT = 1.0 | (0.1, 50.0)    # [m]
  VCmax25::FT = 89.45 | (5.0, 200.0)    # [μmol m-2 s-1]
end
```

---

## 3. 函数定义

### 3.1 文档字符串

```julia
"""
    InitParam_Veg(lc::Int=1; FT=Float64)

读取 JSON 配置文件并返回 ParamVeg 结构体。

# Arguments
- `lc::Int=1`: 植被类型代码 (1=ENF, 6=DBF, 9=EBF)
- `FT::Type=Float64`: 浮点类型

# Returns
- `ParamVeg{FT}`: 植被参数结构体

# Examples
```julia
veg = InitParam_Veg(6)  # DBF
```
"""
function InitParam_Veg(lc::Int=1; FT=Float64)
  # ...
end
```

### 3.2 参数解包

```julia
using UnPack

# 使用 @unpack 提取字段
function some_function(model::ParamBEPS)
  @unpack α_canopy_vis, α_canopy_nir, VCmax25 = model.veg
  @unpack θ_sat, K_sat = model.hydraulic

  # 使用解包的变量
  result = α_canopy_vis * θ_sat
end
```

### 3.3 类型注解

```julia
# 明确类型注解
function photosynthesis_jl(T_leaf_p::T, Rsn_leaf::T, ea::T,
  gb_w::T, Vcmax25::T,
  β_soil::T, g0_w::T, g1_w::T, cii::T,
  T_leaf::T, LH_leaf::T, ca::T=CO2_air) where {T<:Cdouble}
  # ...
end

# 使用类型参数
function InitParam_Soil(SoilType::Int, N::Int, FT::Type)
  hydraulic = ParamSoilHydraulicLayers{FT,N}()
  return hydraulic
end
```

---

## 4. 性能优化

### 4.1 循环优化

```julia
# 使用 @inbounds 避免边界检查
@inbounds for i in 1:n
  θ[i] += (inf - Ett[i]) * Δt / dz[i]
end

# 使用 @fastmath 加速数学运算
@fastmath function cal_K(θ, θ_sat, K_sat, b)
  return K_sat * (θ / θ_sat)^(2b + 3)
end
```

### 4.2 原地操作

```julia
# 使用 .= 进行原地赋值
θ_prev .= θ

# 使用 ! 表示原地修改
function UpdateSoilMoisture!(st::StateBEPS)
  # 直接修改 st
  st.θ .= new_values
end
```

### 4.3 避免不必要的类型转换

```julia
# 不好
for i in 1:n
  value = Float64(array[i])  # 每次都转换
end

# 好
FT = Float64
for i in 1:n
  value = FT(array[i])  # 一次性转换
end
```

---

## 5. 数值安全

### 5.1 边界截断

```julia
# 温度截断
T = clamp(T, -50, 50)
# 或相对于气温
T = clamp(T, Ta-25, Ta+25)

# 土壤水分截断
θ = clamp(θ, θ_vwp, θ_sat)

# 热通量截断
G = clamp(G, -200, 200)

# 冰比例截断
ice_ratio = clamp(ice_ratio, 0.0, 1.0)
```

### 5.2 边界保护

```julia
# 数组访问保护
val_T = i <= length(Tsoil_c) ? Tsoil_c[i] : Tsoil_c[end]

# 条件保护
PPFD = 4.55 * 0.5 * Rsn_leaf
(2PPFD < 1) && (PPFD = 0.0)  # 夜间光合作用为 0
```

### 5.3 物理约束

```julia
# 水分不能为负
θ = max(θ, 0.0)

# 比率在 [0, 1]
fraction = clamp(fraction, 0.0, 1.0)

# 避免除零
denominator = x + eps(x)  # 添加机器精度
```

---

## 6. 代码组织

### 6.1 文件结构

```
src/
├── BEPS.jl                 # 主模块
├── beps_main.jl            # 顶层接口
├── inter_prg.jl            # 主积分循环
├── DataType/
│   ├── Params/
│   │   ├── Params.jl       # 参数结构体
│   │   ├── Param_Init.jl   # 参数初始化
│   │   └── GlobalData.jl   # 全局数据
│   └── DataType.jl         # 类型定义
├── SoilPhysics/
│   ├── UpdateSoilMoisture.jl
│   └── UpdateHeatFlux.jl
└── photosynthesis.jl
```

### 6.2 导入顺序

```julia
# 1. 标准库
using UnPack
using Statistics

# 2. 第三方包
using Parameters
using DataFrames
using FieldMetadata: @metadata

# 3. 本地模块
include("DataType.jl")
include("photosynthesis.jl")
```

### 6.3 函数组织

```julia
# 1. 导出函数
export InitParam_Veg, InitParam_Soil

# 2. 辅助函数
function helper1()
end

function helper2()
end

# 3. 主函数
function main()
  helper1()
  helper2()
end
```

---

## 7. 注释规范

### 7.1 文档字符串（英文）

```julia
"""
Calculate soil hydraulic conductivity using Campbell's equation.

# Arguments
- `θ`: volumetric soil water content [m³/m³]
- `θ_sat`: saturated water content [m³/m³]
- `K_sat`: saturated hydraulic conductivity [m/s]
- `b`: Campbell's b parameter [-]

# Returns
- Hydraulic conductivity [m/s]

# References
- Campbell, 1974
"""
function cal_K(θ, θ_sat, K_sat, b)
  return K_sat * (θ / θ_sat)^(2b + 3)
end
```

### 7.2 行内注释（中文）

```julia
# 计算最大入渗率
# K_sat * (1 + (θ_sat - θ_prev) / dz * ψ_sat / θ_sat * b)
inf_max = f_water[1] * K_sat[1] * (1 + (θ_sat[1] - θ_prev[1]) / dz[1] *
          ψ_sat[1] * b[1] / θ_sat[1])

# 为了解决相互依赖的关系，循环寻找稳态
@inbounds while total_t < kstep
  # ...
end
```

### 7.3 分节注释

```julia
# ===== 1. 参数提取和计算 =====
@unpack α_canopy_vis, α_canopy_nir = ps.veg

# ===== 2. 气象变量初始化 =====
(; Srad, Tair, RH, wind) = forcing

# ===== 3. 表面状态初始化 =====
init_leaf_dbl(Tc_old, Tair - 0.5)

# ===== 4. 亚小时循环 =====
@inbounds for k = 2:kloop+1
  # ...
end
```

---

## 8. 单位标注

### 8.1 注释中标注单位

```julia
# 温度 [°C]
Tair::Float64 = 20.0

# 长度 [m]
z_canopy_o::FT = 1.0 | (0.1, 50.0)

# 通量 [W/m²]
G::Vector{Float64} = zeros(10)

# 速率 [μmol m-2 s-1]
VCmax25::FT = 89.45 | (5.0, 200.0)

# 比率 [-]
f_soilwater::Cdouble = 0
```

### 8.2 物理量计算

```julia
# 明确单位转换
K_sat = FT.(p.K_sat[1:n])    # [m s-1]
K_sat_cm_h = K_sat * 360000   # [cm h-1]

# 能量单位
Lv = 2.501e6  # [J/kg], 汽化潜热
Rn = 500.0    # [W/m²], 净辐射
```

---

## 9. 错误处理

### 9.1 断言

```julia
# 参数验证
@assert 1 <= SoilType <= 11 "SoilType must be between 1 and 11"
@assert 0.0 <= θ <= θ_sat "Soil moisture exceeds saturation"

# 物理约束
@assert value >= 0 "Value must be non-negative"
```

### 9.2 边界检查

```julia
# 使用 clamp 而不是断言
θ = clamp(θ, θ_vwp, θ_sat)

# 友好的错误提示
if !isfinite(value)
  @warn "Non-finite value detected" value=value
  value = 0.0  # 使用默认值
end
```

---

## 10. 测试规范

### 10.1 测试文件命名

```
test/
├── test-beps_main.jl       # 集成测试
├── test-ModelParams.jl     # 参数测试
└── modules/
    ├── test-Soil.jl        # 土壤模块测试
    └── test-photosynthesis.jl
```

### 10.2 测试结构

```julia
using BEPS
using Test

@testset "Soil Physics" begin
  @testset "UpdateSoilMoisture" begin
    soil = init_soil()
    UpdateSoilMoisture(soil, kstep)

    @test all(soil.θ .>= soil.θ_vwp)
    @test all(soil.θ .<= soil.θ_sat)
  end
end
```

---

## 11. 性能基准

### 11.1 使用 BenchmarkTools

```julia
using BenchmarkTools

@btime UpdateSoilMoisture($soil, $kstep)
# 输出: 1.234 μs (100 allocations: 5.67 KiB)
```

### 11.2 性能目标

- 首次运行（含 JIT）: < 5 秒
- 后续运行: < 0.5 秒
- 内存分配: 最小化堆分配

---

## 12. 代码审查清单

提交代码前检查：

- [ ] 函数有文档字符串
- [ ] 物理量有单位标注
- [ ] 数值有边界截断
- [ ] 循环使用 `@inbounds`
- [ ] 类型注解明确
- [ ] 变量命名符合规范
- [ ] 注释清晰（中文行内，英文文档）
- [ ] 测试覆盖核心功能
- [ ] 性能无明显退化

---

## 13. 示例对比

### 不好 vs 好

```julia
# ❌ 不好
function calc(a, b, c)
  x = a * b
  y = x / c
  return y
end

# ✅ 好
"""
Calculate hydraulic conductivity using Campbell's equation.

# Arguments
- `θ`: volumetric water content [m³/m³]
- `θ_sat`: saturated water content [m³/m³]
- `K_sat`: saturated conductivity [m/s]
- `b`: Campbell parameter [-]

# Returns
- Hydraulic conductivity [m/s]
"""
@fastmath function cal_K(θ::T, θ_sat::T, K_sat::T, b::T) where {T<:AbstractFloat}
  return K_sat * (θ / θ_sat)^(2b + 3)
end
```

---

**遵循这些规范，代码将更易读、更安全、更高效。**
