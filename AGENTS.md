# BEPS.jl - AI Agent Guide

## 0. 关键指令

- **语言:** 请使用中文进行对话
- **代码:** 变量名、注释和文档字符串使用英文
- **代码风格:** 遵循 Linux 极简主义；同时遵循代码规范、排版美观

---

## 1. 项目概述

**BEPS.jl** 是 Boreal Ecosystem Productivity Simulator (BEPS，北方生态系统生产力模拟器) 的 Julia 实现版本，是一个陆地表面模型，用于模拟植被光合作用、能量平衡和土壤过程。

- **代码仓库:** https://github.com/jl-pkgs/BEPS.jl
- **版本:** 0.1.8
- **作者:** Dongdong Kong <kongdd@users.noreply.github.com>
- **许可证:** 见 LICENSE.md
- **性能:** 比传统 C 版本快约 2.5 倍
  - Julia: 0.286s (822.38k allocations: 22.998 MiB)
  - C: 0.787s (629.95k allocations: 13.915 MiB)

### 1.1 BEPS 模拟内容

BEPS 以小时为时间步长，模拟以下过程：

- **冠层过程:** 光合作用（区分阳生/阴生叶片）、气孔导度、叶片能量平衡
- **辐射传输:** 太阳辐射截获与分配（直射/散射、可见光/近红外）
- **能量平衡:** 感热、潜热（蒸腾、蒸发）、净辐射
- **土壤过程:** 多层土壤水分（Richards 方程）、土壤温度（热传导）、冻融循环
- **积雪:** 积雪累积、融化和升华（三层方案）
- **水文:** 降水截留、入渗、径流、根系吸水

---

## 2. 技术栈

### 2.1 语言与运行时
- **主要语言:** Julia
- **支持版本:** 1.8, 1.9, 1.10, 1.11
- **操作系统:** Windows（主要开发环境，因为依赖 `BEPS.clang` DLL）

### 2.2 主要依赖（来自 Project.toml）

| 包名 | 版本 | 用途 |
|------|------|------|
| ComponentArrays | 0.15.29 | 参数的结构化数组 |
| DataFrames | 1.5+ | 表格数据处理 |
| DelimitedFiles | 1.9 | 文件读写 |
| DocStringExtensions | 0.9 | 文档支持 |
| FieldMetadata | 0.3 | 结构体字段元数据 |
| Functors | 0.5.2 | 函数映射 |
| JSON | 1.3+ | 配置文件 |
| ModelParams | 0.1.0 | 参数处理 |
| Parameters | 0.12 | 带默认值的结构体 (@with_kw) |
| StaticArrays | 1.9.16 | 快速小型数组 |
| Statistics | 1.11.1 | 统计函数 |
| UnPack | 1.0 | @unpack/@pack! 宏 |
| Libdl | 标准库 | C 接口 DLL 加载 |

### 2.3 外部依赖
- **BEPS.c:** 用于验证的原始 C 实现（可选，位于 `deps/`）
- **libbeps.dll:** 用于对比测试的编译 C 库（仅 Windows）

---

## 3. 项目结构

```
BEPS.jl/
├── Project.toml              # Julia 包配置
├── Manifest.toml             # 依赖锁定文件
├── README.md                 # 用户文档
├── AGENTS.md                 # 本文件 - AI 智能体指南
├── CLAUDE.md                 # 详细技术文档
├── LICENSE.md                # 许可证
│
├── src/                      # 源代码
│   ├── BEPS.jl              # 主模块 - 入口点
│   ├── beps_main.jl         # 顶层接口（besp_main 函数）
│   ├── beps_modern.jl       # 现代 API 接口
│   ├── beps_optimize.jl     # 参数优化工具
│   ├── BEPS_modules.jl      # 模块导出和包含
│   ├── inter_prg.jl         # 主积分循环（核心物理）
│   │
│   ├── DataType/            # 数据结构和类型
│   │   ├── DataType.jl      # 主要类型定义
│   │   ├── BEPS_State.jl    # StateBEPS 和 Soil 结构体
│   │   ├── OUTPUT.jl        # Results 和 OutputET 结构体
│   │   ├── Met.jl           # 气象输入结构体
│   │   ├── Constant.jl      # 物理常数
│   │   ├── CanopyLayer.jl   # 叶片和冠层类型
│   │   ├── setup.jl         # 模型初始化
│   │   └── Params/          # 参数定义
│   │       ├── Params.jl           # ParamVeg、ParamSoil 结构体
│   │       ├── BEPS_Param.jl       # 组合参数结构体
│   │       ├── Param_Init.jl       # 参数初始化
│   │       ├── GlobalData.jl       # 全局植被数据
│   │       └── data/               # JSON 参数数据库
│   │
│   ├── SPAC/                # 土壤-植物-大气连续体
│   │   ├── SPAC.jl          # 太阳几何、气象学
│   │   ├── Leaf.jl          # 叶片尺度计算
│   │   ├── VCmax.jl         # 最大羧化速率
│   │   ├── lai2.jl          # LAI 分配
│   │   └── helper.jl        # 辅助函数
│   │
│   ├── SoilPhysics/         # 土壤物理模块
│   │   ├── SoilPhysics.jl   # 模块导出
│   │   ├── UpdateHeatFlux.jl       # 热传导
│   │   ├── UpdateSoilMoisture.jl   # 水分运移（Richards 方程）
│   │   └── soil_water_factor_v2.jl # 水分胁迫
│   │
│   ├── clang/               # C 接口（仅 Windows）
│   │   ├── BEPS_c.jl        # C 函数包装器
│   │   ├── SOIL_c.jl        # C 土壤结构体
│   │   └── module.jl        # C 模块导出
│   │
│   ├── photosynthesis.jl    # 光合作用模型（Farquhar）
│   ├── photosynthesis_helper.jl
│   ├── netRadiation.jl      # 辐射传输
│   ├── aerodynamic_conductance.jl  # 湍流交换
│   ├── heat_H_and_LE.jl     # 感热/潜热
│   ├── surface_temperature.jl      # 地表能量平衡
│   ├── evaporation_canopy.jl       # 冠层蒸发
│   ├── evaporation_soil.jl         # 土壤蒸发
│   ├── rainfall_stage.jl    # 降水处理
│   ├── snowpack.jl          # 积雪模型（三个阶段）
│   └── backup/              # 废弃代码
│
├── test/                    # 测试套件
│   ├── runtests.jl          # 测试运行器
│   ├── test-beps_main.jl    # 集成测试
│   ├── test-beps_modern.jl  # 现代 API 测试
│   ├── test-ModelParams.jl  # 参数测试
│   └── modules/             # 单元测试
│       ├── modules.jl       # 模块测试运行器
│       ├── test-Soil.jl     # 土壤物理测试
│       ├── test-photosynthesis.jl
│       ├── test-radiation.jl
│       ├── test-snowpack.jl
│       └── ...
│
├── examples/                # 使用示例
│   ├── example_01.qmd       # 基本用法（Quarto）
│   ├── example_02_CUG.qmd   # CUG 站点示例
│   ├── Figure1_compare_with_C.jl
│   ├── Figure3_snowpack_diagnose.jl
│   └── input/               # 示例输入数据
│
├── docs/                    # 文档
│   ├── modules/             # 模块文档（MD）
│   ├── reference/           # 研究论文
│   ├── images/              # 图片
│   └── *.typ                # Typst 源文件
│
├── data/                    # 序列化测试数据
├── deps/                    # 外部依赖（BEPS.c）
└── .github/workflows/       # CI/CD 配置
    └── CI.yml               # GitHub Actions CI
```

---

## 4. 构建与开发命令

### 4.1 环境设置

```bash
# 克隆仓库
git clone https://github.com/jl-pkgs/BEPS.jl
cd BEPS.jl

# 用于与 C 版本对比（仅开发者）
cd deps
git clone https://github.com/jl-pkgs/BEPS.c
```

### 4.2 运行测试

```bash
# 激活项目
julia --project

# 运行所有测试
julia --project -e "using Pkg; Pkg.test()"

# 仅编译检查
julia --project -e "using BEPS"

# 运行特定测试文件
julia --project test/test-beps_main.jl

# 运行模块测试
julia --project test/modules/test-Soil.jl
```

### 4.3 运行模型

```julia
using BEPS

# 加载气象数据
d = deserialize("data/p1_meteo")
lai = readdlm("examples/input/p1_lai.txt")[:]

# 设置参数
par = (
    lon=120.5, lat=30.5,
    VegType=25, SoilType=8,
    clumping=0.85,
    Tsoil0=2.2, θ0=0.4115, z_snow0=0.0
)

# 运行模拟
df_jl, df_ET_jl, Tg, θ = besp_main(d, lai; par..., version="julia")
```

### 4.4 运行示例

```bash
# 运行 Quarto 示例
julia --project examples/example_01.qmd
```

---

## 5. 代码风格与规范

### 5.1 命名规范

| 模式 | 含义 | 示例 |
|------|------|------|
| 希腊字母 | 物理量 | `θ`（水分）、`ψ`（势）、`κ`（导热率）、`ρ`（密度） |
| `_p` 后缀 | 上一时间步 | `Tsoil_p` |
| `_c` 后缀 | 当前时间步 | `Tsoil_c` |
| `_sat` 后缀 | 饱和值 | `θ_sat`、`K_sat` |
| `_vwp` 后缀 | 萎蔫点 | `θ_vwp` |
| `_vfc` 后缀 | 田间持水量 | `θ_vfc` |
| `f_` 前缀 | 分数/因子 | `f_soilwater`、`f_root` |
| `r_` 前缀 | 速率 | `r_rain_g`、`r_drainage` |
| `z_` 前缀 | 深度 | `z_snow`、`z_water` |

**重要:** 使用 `T` 或 `Tsoil` 表示温度。**不要使用 `temp`。**

### 5.2 Julia 编程模式

```julia
# 带默认值的结构体 - 使用 @with_kw
@with_kw mutable struct Soil
    θ::Vector{Float64} = zeros(10)
    Tsoil_c::Vector{Float64} = zeros(10)
end

# 提取字段 - 使用 @unpack
@unpack θ, Tsoil_c, dz = soil

# 赋值回结构体 - 使用 @pack!
@pack! soil = θ, Tsoil_c

# 性能关键循环 - 使用 @inbounds
@inbounds for i = 1:n
    a[i] = b[i] + c[i]
end

# 原地操作 - 使用 .=
soil.θ[1:5] .= new_values

# 适当使用快速数学运算
@fastmath sqrt(x)
```

### 5.3 安全性：数值截断

务必将物理量截断到安全范围内：

```julia
# 温度
clamp(T, -50, 50)
# 或相对于气温
clamp(T, Ta-25, Ta+25)

# 土壤水分
clamp(θ, θ_vwp, θ_sat)

# 热通量
clamp(G, -200, 200)

# 冰比例
clamp(ice_ratio, 0.0, 1.0)
```

### 5.4 文档规范

- **文档字符串:** 使用英文
- **行内注释:** 复杂逻辑可使用/鼓励使用中文
- **参考文献:** 引用论文（如 "Campbell 1974"、"Chen 2007 Eq. 18"）

---

## 6. 架构与数据结构

### 6.1 核心数据流

```
气象输入 (Met)
    ↓
inter_prg_jl() - 主积分循环
    ↓
┌─────────────────┬─────────────────┬─────────────────┐
↓                 ↓                 ↓                 ↓
积雪           辐射             光合作用           土壤物理
(rainfall_     (netRadiation_   (photosynthesis_  (UpdateSoilMoisture、
 stage.jl)     jl)              jl)               UpdateHeatFlux)
    ↓                 ↓                 ↓                 ↓
地表           冠层蒸散发        GPP/蒸腾          土壤温度/水分
状态           (evaporation_)  (transpiration_)  (Root_Water_Uptake)
    ↓
结果 (GPP、ET、SH、LH) + OutputET
```

### 6.2 关键结构体

#### StateBEPS (`src/DataType/BEPS_State.jl`)
跨时间步持续存在的状态变量：

```julia
@with_kw mutable struct StateBEPS
    n_layer::Cint = 5
    dz::Vector{Float64} = zeros(10)      # 层厚度 [m]
    
    # 温度状态
    Tsnow_c::SnowLand{FT}                 # 当前积雪温度
    Tsnow_p::SnowLand{FT}                 # 上一时间步积雪温度
    Tsoil_c::Vector{Float64} = zeros(10)  # 当前土壤温度 [°C]
    Tsoil_p::Vector{Float64} = zeros(10)  # 上一时间步土壤温度 [°C]
    
    # 水分状态
    θ::Vector{Float64} = zeros(10)        # 土壤水分 [m³/m³]
    θ_prev::Vector{Float64} = zeros(10)   # 上一时间步水分
    z_water::Cdouble = 0                  # 积水深度 [m]
    z_snow::Cdouble = 0                   # 积雪深度 [m]
    
    # 物理属性（每步更新）
    κ::Vector{Float64} = zeros(10)        # 热导率 [W/m/K]
    Cv::Vector{Float64} = zeros(10)       # 热容量 [J/m³/K]
    ψ::Vector{Float64} = zeros(10)        # 基质势 [m]
    ice_ratio::Vector{Float64} = zeros(10)# 冻结比例 [0-1]
    
    # 通量
    G::Vector{Float64} = zeros(10)        # 热通量 [W/m²]
    Ett::Vector{Float64} = zeros(10)      # 每层蒸散发 [m/s]
    
    # 水分胁迫
    f_soilwater::Cdouble = 0              # 整体水分胁迫因子
    f_root::Vector{Float64} = zeros(10)   # 每层根系比例
end
```

#### ParamBEPS (`src/DataType/Params/BEPS_Param.jl`)
模拟期间保持不变的模型参数：

```julia
@with_kw mutable struct ParamBEPS{FT}
    veg::ParamVeg{FT}                     # 植被参数
    hydraulic::ParamSoilHydraulicLayers{FT}  # 分层的土壤水力参数
    thermal::ParamSoilThermalLayers{FT}      # 分层的土壤热参数
end
```

#### Met (`src/DataType/Met.jl`)
气象强迫：

```julia
@with_kw mutable struct Met
    Srad::Cdouble = 0    # 短波辐射 [W/m²]
    Lrad::Cdouble = 0    # 长波辐射 [W/m²]
    Tair::Cdouble = 0    # 气温 [°C]
    hum::Cdouble = 0     # 相对湿度 [%]
    rain::Cdouble = 0    # 降水 [mm/hour]
    wind::Cdouble = 0    # 风速 [m/s]
end
```

### 6.3 积分循环 (`src/inter_prg.jl`)

主模拟循环在每个次小时时间步（360秒，每小时10步）执行以下步骤：

```julia
for k = 2:kloop+1  # 10 个次小时步
    # 1. 积雪阶段 1 - 降水截留
    z_snow = snowpack_stage1_jl(...)
    
    # 2. 降雨阶段 1 - 冠层截留
    r_rain_g = rainfall_stage1_jl(...)
    
    # 3. 土壤水分胁迫因子
    soil_water_factor_v2(state, ps)
    
    # 4. 冠层能量平衡迭代
    radiation_o, radiation_u, radiation_g, ra_g = 
        solve_canopy_energy_balance!(...)
    
    # 5. 蒸腾和蒸发
    Trans_o, Trans_u = transpiration_jl(...)
    Eil_o, Eil_u, EiS_o, EiS_u = evaporation_canopy_jl(...)
    
    # 6. 降雨阶段 2 - 更新冠层水分
    rainfall_stage2_jl(...)
    
    # 7. 积雪阶段 2 - 更新冠层积雪
    snowpack_stage2_jl(...)
    
    # 8. 土壤蒸发
    Evap_soil, Evap_SW, Evap_SS = evaporation_soil_jl(...)
    
    # 9. 地表温度（能量平衡耦合）
    G[1] = surface_temperature!(...)
    
    # 10. 积雪阶段 3 - 融雪/冻结
    z_snow, z_water = snowpack_stage3_jl(...)
    
    # 11. 感热通量
    Qhc_o, Qhc_u, Qhg = sensible_heat_jl(...)
    
    # 12. 土壤热通量（第2-5层）
    UpdateHeatFlux(state, Tair, kstep)
    
    # 13. 根系吸水分配
    Root_Water_Uptake(state, Trans_o, Trans_u, Evap_soil)
    
    # 14. 土壤水分更新（Richards 方程）
    UpdateSoilMoisture(state, ps, kstep)
end
```

### 6.4 常数 (`src/DataType/Constant.jl`)

关键物理常数：
- `step = 3600.0` - 小时时间步 [s]
- `kstep = 360.0` - 次小时时间步 [s]
- `kloop = 10` - 每小时次小时迭代次数
- `layer = 5` - 土壤层数
- `CO2_air = 380.0` - 大气 CO₂ [ppm]
- `ρₐ = 1.292` - 空气密度 [kg/m³]
- `ρ_w = 1025.0` - 水密度 [kg/m³]
- `Lv_solid = 2.83e6` - 升华潜热 [J/kg]

---

## 7. 测试策略

### 7.1 测试框架
- **框架:** Julia 内置 `Test` 模块
- **验证:** 主要通过与 C 实现（`BEPS.c`）对比验证
- **结构:** 单元测试 + 集成测试

### 7.2 测试组织

```
test/
├── runtests.jl              # 主测试运行器
├── test-beps_main.jl        # 集成测试（Julia vs C）
├── test-beps_modern.jl      # 现代 API 测试
├── test-ModelParams.jl      # 参数处理测试
├── test-utilize.jl          # 工具函数测试
└── modules/                 # 单元测试
    ├── test-Soil.jl         # 土壤物理（详细 C 对比）
    ├── test-photosynthesis.jl
    ├── test-radiation.jl
    ├── test-snowpack.jl
    ├── test-sensible_heat.jl
    ├── test-aerodynamic_conductance.jl
    ├── test-rainfall_stage1.jl
    ├── test-param.jl
    └── test-setup.jl
```

### 7.3 关键测试模式

```julia
# 对比 Julia 与 C 实现
@testset "UpdateHeatFlux" begin
    p_jl, p_c = init_soil()
    
    # 运行两种实现
    UpdateHeatFlux(p_jl, Ta, kstep)
    UpdateHeatFlux(p_c, Ta, kstep)
    
    # 对比所有字段
    is_soil_equal(p_jl, p_c; tol=1e-7)
end
```

### 7.4 验证检查

1. **质量平衡:** `Σ(Δθ * dz) = 入渗 - 径流 - ET - 排水`
2. **能量平衡:** `Rn - LE - H - G ≈ 0`
3. **物理边界:**
   - `θ ∈ [θ_vwp, θ_sat]`
   - `T ∈ [-50, 50]` °C
   - `ice_ratio ∈ [0, 1]`
4. **C/Julia 对比:** 所有土壤结构体字段在容差内匹配

### 7.5 运行测试

```bash
# 所有测试
julia --project -e "using Pkg; Pkg.test()"

# 特定测试集
julia --project -e "using Pkg; Pkg.test(test_args=\"Soil\")"

# 单独测试文件
julia --project test/modules/test-Soil.jl
```

---

## 8. CI/CD 配置

### 8.1 GitHub Actions (`.github/workflows/CI.yml`)

```yaml
- 触发条件: push 到 src/** 或 test/**，pull requests
- Julia 版本: 1（最新稳定版）
- 操作系统: windows-latest（主要）
- 步骤:
  1. 检出代码
  2. 设置 Julia
  3. 构建包
  4. 运行测试
  5. 处理覆盖率
  6. 上传到 Codecov
```

### 8.2 覆盖率
- **服务:** Codecov
- **配置:** `codecov.yml`
- **目标:** 保持土壤物理模块的高覆盖率

---

## 9. 文档

### 9.1 文档文件

| 文件 | 内容 |
|------|------|
| `README.md` | 面向用户的快速入门和概览 |
| `AGENTS.md` | 本文件 - AI 编程智能体指南 |
| `CLAUDE.md` | 详细技术文档 |
| `docs/modules/*.md` | 模块特定文档 |
| `docs/reference/` | 研究论文和参考文献 |

### 9.2 关键参考文献

- **Campbell (1974)** - 土壤水分保持和水力传导
- **Chen (2007)** - 体积热容量（公式 18）
- **Bonan (2019)** - Campbell 参数表 8.2
- **He et al. (2017)** - 土壤水分胁迫因子（JGR-B，公式 4-5）
- **Farquhar et al.** - 光合作用模型
- **BEPS V2023** - 当前实现版本

---

## 10. 开发工作流

### 10.1 进行修改

1. **理解物理:** 查阅 CLAUDE.md 和 docs/modules/
2. **保持 C 兼容性:** 如果修改核心物理，更新 C 包装器
3. **添加测试:** 确保通过与 C 或解析解对比验证变更
4. **更新文档:** 保持 CLAUDE.md 和 AGENTS.md 同步
5. **运行完整测试套件:** 提交前

### 10.2 Bug 修复流程

参见 README.md "Bugs Fixed" 部分了解示例：

```markdown
### YYYY-MM-DD
- [x] `function_name`: Bug 描述
  - 根本原因
  - 应用的修复
  - 添加的测试
```

### 10.3 添加新功能

1. 在 `DataType/` 中使用 `@with_kw` 定义新结构体
2. 在 `setup.jl` 或 `Param_Init.jl` 中添加初始化
3. 在相应模块中实现核心逻辑
4. 在 `inter_prg.jl` 中添加至积分循环
5. 在 `BEPS_modules.jl` 中导出
6. 在 `test/modules/` 中编写测试
7. 在 `docs/modules/` 中编写文档

---

## 11. 已知问题与限制

### 11.1 平台支持
- **Windows:** 完全支持（包括 C 对比）
- **Linux/macOS:** 仅 Julia 模式（无 C DLL）

### 11.2 待解决问题
- 植被光合作用的叶片温度参数传递需要审查

### 11.3 性能说明
- 首次运行包含 JIT 编译开销
- 使用 `--project` 保持环境一致
- 基准测试前使用 `using BEPS` 预编译

---

## 12. 快速参考

### 12.1 初始化模型

```julia
using BEPS

# 设置
soil, state, params = setup_model(VegType, SoilType;
    Ta, Tsoil, θ0, z_snow, r_drainage, r_root_decay)
```

### 12.2 运行单个时间步

```julia
# 1. 土壤水分胁迫
soil_water_factor_v2(state, ps)

# 2. 冠层过程（辐射、光合作用、ET）
# ... 见 inter_prg.jl

# 3. 热属性
UpdateThermal_κ(state)
UpdateThermal_Cv(state)

# 4. 地表温度
G = surface_temperature!(state, ps, prev, curr, ...)

# 5. 热通量
UpdateHeatFlux(state, Ta, kstep)

# 6. 根系吸水
Root_Water_Uptake(state, Trans_o, Trans_u, Evap_soil)

# 7. 水分更新
UpdateSoilMoisture(state, ps, kstep)
```

### 12.3 访问结果

```julia
θ = state.θ[1:5]           # 分层土壤水分
T = state.Tsoil_c[1:5]     # 分层土壤温度
ice = state.ice_ratio[1:5] # 分层冰比例
G = state.G[1:5]           # 分层热通量
```

---

## 13. 联系与资源

- **代码仓库:** https://github.com/jl-pkgs/BEPS.jl
- **文档:** https://jl-pkgs.github.io/BEPS.jl/dev
- **问题:** https://github.com/jl-pkgs/BEPS.jl/issues

**记住:** 如有疑问，请参考 `deps/BEPS.c` 中的 C 实现作为权威物理依据。
