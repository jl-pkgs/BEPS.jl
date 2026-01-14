# 地表温度模块文档

本文档详细描述了 `src/surface_temperature.jl` 中的 `surface_temperature_jl` 函数及其相关的数据结构、物理模型和边界条件。

## 1. 核心函数概览

### `surface_temperature_jl`

该函数基于能量平衡原理，计算地表各层温度（混合地表、土壤表面、积雪层）以及向下传入土壤的热通量。

```julia
function surface_temperature_jl(T_air, rh_air, z_snow, z_water,
  cp_soil1, cp_soil0, Gheat_g,
  z_soil1, ρ_snow, tempL_u, Rn_g,
  E_soil, E_water_g, E_snow_g, λ_soil1,
  perc_snow_g, G_soil1,
  T_ground_last,
  T_soil1_last, T_soil_surf_last, T_soil0_last,
  T_snow_last, T_snow1_last, T_snow2_last)
```

**主要输入/输出:**

*   **输入**: 气象数据 (`T_air`, `rh_air`), 积雪状态 (`z_snow`, `ρ_snow`), 土壤热参数 (`cp`, `λ`), 辐射与蒸散 (`Rn_g`, `E_...`), 上一时刻状态 (`_last`).
*   **输出**: 
    *   `G`: 地表热通量 [W m-2]
    *   `T_ground`: 混合地表温度 [°C]
    *   `T_soil_surf`: 土壤表面温度 [°C]
    *   `T_soil0`: 有效表层土壤温度 [°C]
    *   `T_snow...`: 积雪各层温度 [°C]

---

## 2. 物理模型与变量定义

### 2.1 关键温度变量辨析

模型区分了三种不同的“地表温度”概念，以处理积雪覆盖的不均匀性。

| 变量名 | 物理含义 | 说明 |
| :--- | :--- | :--- |
| **`T_ground`** | **混合地表温度** (Bulk Surface Temp) | 用于计算向上的长波辐射。代表整个网格单元（雪面+裸土）的平均辐射温度。 |
| **`T_soil_surf`** | **平均土壤表面温度** (Avg Soil Surface Temp) | 土壤层的上边界温度。如果是部分积雪，它是裸土表面温度和雪下土壤表面温度的**加权平均值**。 |
| **`T_soil0`** | **裸土/有效表层温度** (Effective Surface Soil Temp) | 特指直接暴露在空气中的裸土表面温度（在斑块雪情况下）。全覆盖时等于 `T_soil_surf`。 |

### 2.2 积雪分层策略

根据积雪深度 `z_snow`，模型采用不同的分层策略：

1.  **无雪/极薄雪 (`z_snow <= 0.02m`)**:
    *   单层模型。
    *   `T_ground` = `T_soil_surf` = `T_soil0`。
    *   直接计算土壤表面能量平衡。

2.  **薄雪/斑块雪 (`0.02m < z_snow <= 0.05m`)**:
    *   双板块模型 (Two-tile approach)。
    *   **裸土部分**: 独立计算 `T_soil0`。
    *   **积雪部分**: 单层积雪模型计算 `T_snow`。
    *   **聚合**: `T_soil_surf` 和 `T_ground` 由两者加权平均得出。

3.  **深雪 (`z_snow > 0.05m`)**:
    *   多层积雪模型 (最多3层节点: Surface, Middle, Bottom)。
    *   积雪层完全隔绝了土壤与大气。
    *   `T_soil_surf` 仅受积雪底层热通量 (`G_snow2`) 和土壤深层热传导 (`G_soil1`) 控制。

## 3. 数学原理与公式推导

### 3.1 通用隐式求解 (Implicit Solver)

对于地表温度或积雪温度，我们通常求解以下形式的能量平衡方程：

$$ C \frac{\partial T}{\partial t} = R_n - H - LE - G $$

离散化为隐式格式：

$$ C \Delta z \frac{T^{n+1} - T^n}{\Delta t} = R_n - \rho C_p \frac{T^{n+1} - T_{air}}{r_a} - \lambda \frac{T^{n+1} - T_{bot}}{z} $$

整理为线性方程 $T^{n+1} = \frac{Numerator}{Denominator}$，即 `solve_imp` 函数：

**公式:**
$$ T_{new} = \frac{T_{old} \cdot I + G_{net} \cdot r_a \cdot z_{rad} + \rho C_p \cdot T_{bnd} \cdot z + c_s \cdot r_a \cdot \lambda_{bot} \cdot T_{bot}}{\rho C_p \cdot z + c_s \cdot r_a \cdot \lambda_{bot} + I} $$

其中：
*   $I = \frac{C \Delta z}{\Delta t} \cdot r_a \cdot z$: 热惯性项
*   $G_{net}$: 净辐射及其他源项
*   $T_{bnd}$: 上边界温度 (通常为 $T_{air}$)
*   $T_{bot}$: 下边界温度 (下一层土壤或雪层)
*   $\lambda_{bot}$: 下层热导率
*   $r_a$: 空气动力学阻力
*   $c_s$: 缩放系数 (用于调整热传导项)

### 3.2 显式时间步进 (Explicit Time Stepping)

对于内部节点（如积雪中间层），使用简单的显式欧拉法 `step_exp`：

$$ C \Delta z \frac{T^{n+1} - T^n}{\Delta t} = F_{in} - F_{out} $$

**公式:**
$$ T_{new} = T_{old} + \frac{F_{in} - F_{out}}{C \Delta z} \cdot \Delta t $$

### 3.3 热导率计算

积雪热导率随密度变化，采用 Jordan (1991) 经验公式 `cal_κ_snow`：

$$ \kappa_{snow} = 0.021 + 4.2 \times 10^{-4} \cdot \rho + 2.2 \times 10^{-9} \cdot \rho^3 $$

---

## 4. 辅助函数与数据结构说明

### 4.1 数据结构维度 (`layer + 2`)

在 `TransientCache` 结构体中，`Cs` (土壤体积热容)、`G` (土壤热通量) 和 `T_soil` (土壤温度) 的维度被定义为 `layer + 2`。

假设实际物理土壤层数为 `N`:

| 索引 (Index) | 对应物理层级 | 变量含义 | 备注 |
| :--- | :--- | :--- | :--- |
| **1** | **Layer 0 (Surface Skin)** | $T_{surf}$, $G_{ground}$ | **地表界面层**。<br>`T`对应 `T_soil_surf`。<br>`G[1]` 是从地表向下传入第一层土壤的通量。 |
| **2** | **Layer 1 (First Soil)** | $T_{soil,1}$, $G_{1 \to 2}$ | **第一层实际土壤**。<br>`T`对应第一层中心温度。<br>`G[2]` 是第一层底部的通量。 |
| **3 ... N+1** | **Layer 2 ... N** | $T_{soil,i}$, $G_{i \to i+1}$ | 后续物理土壤层。 |
| **N+2** | **Bottom Boundary** | $G_{bottom}$ | **下边界**。<br>存储最底层下方的边界通量。 |

### 4.2 索引偏移说明

在辅助函数 `surface_temperature!` 中：
*   `κ[2]` 被赋值为 `soil.κ[1]`，代表第一层土壤的热导率。
*   调用 `surface_temperature_jl` 时，`κ[2]` 作为 `λ_soil1` 传入。

### 4.3 边界条件

*   **上边界**: `G[1]` 由 `surface_temperature_jl` 计算，作为源项驱动土壤热扩散。
*   **下边界**: `G[N+2]` 采用恒温边界条件，即假设深层土壤温度恒定为年平均气温。
