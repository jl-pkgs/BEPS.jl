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

## 实现细节与参数映射

### 索引偏移说明

在 `surface_temperature_helpers.jl` 的 `prepare_and_compute_surface_temperature!` 函数中，存在一个索引偏移，需要注意：

*   **`κ[2]` 对应 `soil.κ[1]`**：
    *   `soil.κ` 存储实际土壤分层的热导率，其中索引 1 代表最上层土壤（Layer 1）。
    *   在辅助函数的临时数组 `κ` 中，索引 1 通常预留给地表“皮肤”层（Layer 0）或边界条件。
    *   因此，`κ[2]` 被赋值为 `soil.κ[1]`，并在调用核心函数 `surface_temperature_jl` 时作为参数 `λ_soil1`（第一层土壤热导率）传入。

代码示例：
```julia
# src/surface_temperature_helpers.jl

κ[2] = soil.κ[1]                  # 热导率 [W m-1 K-1]
# ...
surface_temperature_jl(...,
  κ[2], # 对应参数 λ_soil1
  ...
)
```

### 数据结构维度说明 (`Cs`, `G`, `T_soil`)

在 `TransientCache` 结构体中，`Cs` (土壤体积热容)、`G` (土壤热通量) 和 `T_soil` (土壤温度) 的维度被定义为 `layer + 2`。这是为了适应多层土壤热传输模型的数值计算需求，并与 Fortran 原始代码的索引习惯保持一致（或兼容）。

假设实际物理土壤层数为 `layer` (通常为 5 或 10)，则数组各层的物理含义如下：

*   **索引 1 (Layer 0 / Surface Skin)**:
    *   **含义**: 地表“皮肤”层（Skin Layer）或界面层。
    *   **T_soil[1]**: 对应 `T_soil_surf` 或 `T_ground`，即与大气直接交换能量的地表温度。
    *   **G[1]**: **地表热通量** ($G_{ground}$)，即从地表传入第一层土壤的热通量。
    *   **Cs[1]**: 地表层的热容（通常取第一层土壤的值或特定参数）。

*   **索引 2 (Layer 1 / First Soil Layer)**:
    *   **含义**: **第一层实际物理土壤**。
    *   **T_soil[2]**: 第一层土壤中心的温度。
    *   **G[2]**: 第一层土壤底部的热通量（或第一层与第二层之间的通量）。
    *   **Cs[2]**: 第一层土壤的体积热容。
    *   *注：在 `surface_temperature_helpers.jl` 中，`κ[2]` 对应 `soil.κ[1]`，也是指向这一层。*

*   **索引 3 ... layer+1 (Layer 2 ... Last Layer)**:
    *   **含义**: 后续的物理土壤层。
    *   **T_soil[i]**: 第 `i-1` 层实际土壤的温度。
    *   **G[i]**: 第 `i-1` 层土壤底部的热通量。

*   **索引 layer+2 (Bottom Boundary)**:
    *   **含义**: 底部边界或计算缓冲区。
    *   **用途**: 用于存储最底层下方的热通量边界条件（如零通量或恒温边界），或者作为循环计算时的 Ghost Point 防止数组越界。

**总结表**:

| 索引 (Index) | 物理层级 | 变量含义 | 备注 |
| :--- | :--- | :--- | :--- |
| **1** | Layer 0 (Surface) | $T_{surf}$, $G_{ground}$ | 地表界面 |
| **2** | Layer 1 (Soil) | $T_{soil,1}$, $G_{1 \to 2}$ | 第一层实际土壤 |
| **3** | Layer 2 (Soil) | $T_{soil,2}$, $G_{2 \to 3}$ | 第二层实际土壤 |
| ... | ... | ... | ... |
| **layer+1** | Layer N (Soil) | $T_{soil,N}$, $G_{N \to bottom}$ | 最后一层实际土壤 |
| **layer+2** | Bottom | Boundary Condition | 底部边界/缓冲 |

### 边界条件处理

#### 上边界 (Upper Boundary)
上边界（地表）的热通量 `G[1]` 由 `surface_temperature_jl` 模块直接计算：
*   **计算依据**: 地表能量平衡（Net Radiation - Sensible Heat - Latent Heat）。
*   **作用**: 作为驱动力输入到土壤热传输方程中，更新第一层土壤温度。

#### 下边界 (Lower Boundary)
下边界处理位于 `src/Soil/UpdateHeatFlux.jl` 中，采用**恒温边界条件**：
*   **假设**: 在最后一层土壤下方一定深度处，温度恒定为年平均气温。
*   **实现**: 计算最后一层土壤与恒温边界之间的热传导通量。
    ```julia
    G[n+1] = κ[n] * (Tsoil_p[n] - Tair_annual_mean) / (DEPTH_F + dz[n] * 0.5)
    ```
    其中 `DEPTH_F` (默认 6m) 是从最后一层中心到恒温边界的额外距离。

```
  大气 (T_air)
    |
    | ra_g (空气动力学阻力)
    |
------------------  <-- T_ground (混合地表温度)
|     积雪层     |      (如果有雪, T_ground ≈ T_snow)
|  [ T_snow1 ]   |  <-- T_snow1 (积雪上层温度)
|----------------|
|  [ T_snow2 ]   |  <-- T_snow2 (积雪下层温度)
------------------  <-- T_soil_surf (土壤表面温度 / 雪-土界面)
    |                 (如果无雪, T_soil_surf ≈ T_ground)
    |
0层土壤           <-- T_soil0 (有效表层土壤温度)
    |
------------------
    |
1层土壤           <-- T_soil1 (第一层土壤温度)
    |
```

## 3. 数据结构与维度细节

### 3.1 数组维度设计 (`layer + 2`)

在 `TransientCache` 结构体中，`Cs`, `G`, `T_soil` 等数组的维度被定义为 `layer + 2`。这是为了完整描述从地表到深层边界的热传输过程。

假设实际物理土壤层数为 `N` (例如 `p.n_layer`):

| 索引 (Index) | 对应物理层级 | 变量含义 | 备注 |
| :--- | :--- | :--- | :--- |
| **1** | **Layer 0 (Surface Skin)** | $T_{surf}$, $G_{ground}$ | **地表界面层**。<br>`T`对应 `T_soil_surf`。<br>`G[1]` 是从地表向下传入第一层土壤的通量。 |
| **2** | **Layer 1 (First Soil)** | $T_{soil,1}$, $G_{1 \to 2}$ | **第一层实际土壤**。<br>`T`对应第一层中心温度。<br>`G[2]` 是第一层底部的通量。 |
| **3 ... N+1** | **Layer 2 ... N** | $T_{soil,i}$, $G_{i \to i+1}$ | 后续物理土壤层。 |
| **N+2** | **Bottom Boundary** | $G_{bottom}$ | **下边界**。<br>存储最底层下方的边界通量。 |

### 3.2 索引偏移说明

在辅助函数 `surface_temperature_helpers.jl` 中存在索引映射：

*   **`soil.κ` (参数)**: 索引 1 代表第一层实际土壤。
*   **`κ` (临时数组)**: 为了匹配上述 `layer+2` 结构，将 `soil.κ[1]` 赋值给 `κ[2]`。
*   **调用**: `surface_temperature_jl` 接收 `κ[2]` 作为 `λ_soil1` 参数。

---

## 4. 边界条件处理

### 4.1 上边界 (Upper Boundary)
*   **位置**: `G[1]` (索引 1)。
*   **驱动**: **地表能量平衡**。
*   **实现**: 由 `surface_temperature_jl` 函数计算。综合了净辐射 (`Rn`)、感热 (`H`)、潜热 (`LE`) 以及积雪热传导。
*   **作用**: 作为源项输入到土壤热扩散方程中，驱动土壤温度变化。

### 4.2 下边界 (Lower Boundary)
*   **位置**: `G[N+2]` (索引 layer+2)。
*   **驱动**: **恒温边界条件 (Dirichlet)**。
*   **实现**: 在 `src/Soil/UpdateHeatFlux.jl` 中计算。
*   **公式**: 
    $$ G_{bottom} = \kappa_{N} \frac{T_{soil,N} - T_{annual\_mean}}{D_F + 0.5 \Delta z_N} $$
    其中 $D_F$ (默认 6m) 是从最后一层底部到恒温层的距离。
