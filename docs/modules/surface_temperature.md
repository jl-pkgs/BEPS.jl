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

### 3.1 隐式求解 (Implicit Solver)

对于地表温度或积雪温度，我们通常求解以下形式的能量平衡方程（针对厚度为 $\Delta z$ 的层）：

$$
C_v \Delta z \frac{\partial T}{\partial t} = G_{net} - H - G_{down} 
$$

其中 $G_{net} = R_n - LE$ 为净辐射及其他源项，$G_{down}$ 为向下传导热通量。

离散化为隐式格式（Implicit Euler）：

$$
\Delta E (T^{n+1} - T^n) = G_{net} \cdot \frac{z_{rad}}{z} - \rho C_p \frac{T^{n+1} - T_{up}}{r_a} - \kappa_{bot} \frac{c_s}{z} (T^{n+1} - T_{bot})
$$

其中 $\Delta E = \frac{C_v \Delta z}{\Delta t}$。引入比例因子 $z_{rad}/z$ 和 $c_s$ 是为了处理不同物理过程的作用深度差异（如辐射吸收深度与热传导距离不同）。

为消除分母，方程两边同乘以 $r_a \cdot z$：

$$
\Delta E \cdot r_a z (T^{n+1} - T^n) = G_{net} r_a z_{rad} - \rho C_p z (T^{n+1} - T_{up}) - c_s r_a \kappa_{bot} (T^{n+1} - T_{bot})
$$

令惯性项 $I = \Delta E \cdot r_a \cdot z$，展开并整理 $T^{n+1}$ 项：

$$
I T^{n+1} - I T^n = G_{net} r_a z_{rad} - \rho C_p z T^{n+1} + \rho C_p z T_{up} - c_s r_a \kappa_{bot} T^{n+1} + c_s r_a \kappa_{bot} T_{bot}
$$

移项合并 $T^{n+1}$：

$$
T^{n+1} (I + \rho C_p z + c_s r_a \kappa_{bot}) = I T^n + G_{net} r_a z_{rad} + \rho C_p z T_{up} + c_s r_a \kappa_{bot} T_{bot}
$$

最终得到线性方程 $T^{n+1} = \frac{Numerator}{Denominator}$，即 `solve_imp` 函数：

**公式:**
$$
T_{new} = \frac{T_{old} \cdot I + G_{net} \cdot r_a \cdot z_{rad} + \rho C_p \cdot T_{up} \cdot z + c_s \cdot r_a \cdot \kappa_{bot} \cdot T_{bot}}{\rho C_p \cdot z + c_s \cdot r_a \cdot \kappa_{bot} + I} 
$$

其中：
*   $I = \frac{C_v \Delta z}{\Delta t} \cdot r_a \cdot z$: 热惯性项 [$J \cdot m^{-2} \cdot K^{-1}$]
*   $C_v$: 介质（土壤或积雪）的**体积热容** [$J \cdot m^{-3} \cdot K^{-1}$]，满足 $C_v = \rho_{medium} C_{p,medium}$
*   $\rho C_p$: **空气**的体积热容 [$J \cdot m^{-3} \cdot K^{-1}$]，用于感热通量计算
*   $G_{net}$: 净辐射及其他源项 [$W \cdot m^{-2}$]
*   $T_{up}$: 上边界温度 (通常为 $T_{air}$ 或 $T_{c,u}$) [$^\circ C$]
*   $T_{bot}$: 下边界温度 (下一层土壤或雪层) [$^\circ C$]
*   $\kappa_{bot}$: 下层热导率 [$W \cdot m^{-1} \cdot K^{-1}$]
*   $r_a$: 空气动力学阻力 [$s \cdot m^{-1}$]
*   $z, z_{rad}, \Delta z$: 厚度或特征距离 [$m$]
*   $c_s$: 缩放系数 [无量纲]

### 3.2 显式求解 (Explicit Time Stepping)

对于内部节点（如积雪中间层），使用简单的显式欧拉法 `step_exp`：

$$
C_v \Delta z \frac{T^{n+1} - T^n}{\Delta t} = F_{in} - F_{out} 
$$

**公式:**
$$
T_{new} = T_{old} + \frac{F_{in} - F_{out}}{C_v \Delta z} \cdot \Delta t 
$$

### 3.3 热导率计算

积雪热导率随密度变化，采用 Jordan (1991) 经验公式 `cal_κ_snow`：

$$ 
\kappa_{snow} = 0.021 + 4.2 \times 10^{-4} \cdot \rho + 2.2 \times 10^{-9} \cdot \rho^3 
$$

---

## 4. 网格系统与边界层详解 (1:N+2)

BEPS 模型在垂直方向上采用节点中心网格系统，为了方便处理边界通量，数组维度被扩展为 `N + 2`，其中 `N` 为实际物理土壤层数。

### 4.1 逐层详细定义

假设实际物理土壤层数为 `N` (例如 `p.n_layer = 10`)，数组索引 `i` 从 `1` 到 `N+2`：

#### **索引 1: 上边界层 (The Upper Boundary / Surface Skin)**
*   **物理含义**: 地表“皮肤”层，即土壤与大气接触的无限薄界面。
*   **温度 $T[1]$**: 代表地表温度 ($T_{skin}$ 或 $T_{surf}$)。在有雪的情况下，这可能是积雪表面温度或雪-土混合温度。
*   **热通量 $G[1]$**: 代表**进入土壤系统的净热通量** ($G_{ground}$)。
    *   计算来源: 地表能量平衡方程 ($R_n - H - LE$)。
    *   作用: 作为第一层物理土壤的热源边界条件 (Neumann 边界类型)。

#### **索引 2 到 N+1: 物理土壤层 (Physical Soil Domain)**
这是实际进行水分和热量存储、传输的区域。
*   **映射关系**: 数组索引 `i` 对应物理层 `i-1`。
    *   索引 `2` $\rightarrow$ **第 1 层土壤** (最表层土壤)。
    *   索引 `3` $\rightarrow$ **第 2 层土壤**。
    *   ...
    *   索引 `N+1` $\rightarrow$ **第 N 层土壤** (最底层土壤)。
*   **温度 $T[i]$**: 第 `i-1` 层土壤节点的中心温度。
*   **热通量 $G[i]$**: 第 `i-1` 层底部界面流出的热通量（即流向下一层的通量）。

#### **索引 N+2: 下边界层 (The Lower Boundary)**
*   **物理含义**: 位于最后一层物理土壤下方的虚拟边界节点。
*   **温度 $T[N+2]$**: 通常用于存储深层恒温值 ($T_{deep} \approx T_{air, annual\_mean}$)。
*   **热通量 $G[N+2]$**: 代表从模型底部边界流失（或流入）的热通量。
    *   **恒温边界 (Dirichlet)**: 假设深处温度恒定，计算 $G[N+2] = \kappa \frac{T[N+1] - T_{deep}}{\Delta z}$。
    *   **零通量边界 (Neumann)**: 假设绝热，$G[N+2] = 0$。
    *   **BEPS 实现**: 采用恒温边界条件，连接到年平均气温。

### 4.2 总结表

| 数组索引 `i` | 物理层级名称 | 温度变量 $T[i]$ | 热通量变量 $G[i]$ | 备注 |
| :--- | :--- | :--- | :--- | :--- |
| **1** | **上边界 (Surface)** | $T_{surf}$ (地表温度) | $G_{in}$ (入渗通量) | 能量平衡计算得出的边界驱动力 |
| **2** | **Soil Layer 1** | $T_{soil,1}$ | $Flux_{1 \to 2}$ | 实际土壤第一层 |
| **3** | **Soil Layer 2** | $T_{soil,2}$ | $Flux_{2 \to 3}$ | 实际土壤第二层 |
| ... | ... | ... | ... | ... |
| **k** | **Soil Layer k-1** | $T_{soil,k-1}$ | $Flux_{k-1 \to k}$ | 中间层 |
| ... | ... | ... | ... | ... |
| **N+1** | **Soil Layer N** | $T_{soil,N}$ | $Flux_{N \to Bottom}$ | 实际土壤最后一层 |
| **N+2** | **下边界 (Bottom)** | $T_{deep}$ (虚拟) | $Flux_{lost}$ (边界通量) | 虚拟层，用于闭合方程组 |

### 4.3 索引偏移注意

在辅助函数 `surface_temperature!` 中存在索引映射：
*   **`soil.κ` (参数)**: 索引 1 代表第一层实际土壤。
*   **`κ` (临时数组)**: 为了匹配 `layer+2` 结构，将 `soil.κ[1]` 赋值给 `κ[2]`。
*   **调用**: `surface_temperature_jl` 接收 `κ[2]` 作为 `λ_soil1` 参数。
