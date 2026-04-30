# Case 3: Deep Snow Surface Temperature Derivation

本文档针对 `src/surface_temperature.jl` 中 **Case 3 (深雪, $z_{snow} > 0.05m$)** 的物理模型进行公式推导与代码修正建议。

## 1. 物理模型与节点定义

**模型结构 (从上至下):**
1.  **Snow L0 (Top)**: $dz_{s1}=0.02m$, $T_{snow}$ (Implicit)
2.  **Snow L1 (Mid)**: $dz_{s2}=0.02m$, $T_{snow1}$ (Explicit)
3.  **Snow L2 (Bot)**: $dz_{s3}=z_{snow}-0.04m$, $T_{snow2}$ (Explicit)
4.  **Soil Surface**: $dz_{soil0}=0.02m$, $T_{soil\_surf}$ (Explicit)

**假设**: 节点位于各层几何中心 (Node-centered scheme)。

## 2. 热传导距离推导

| 路径 | 节点深度计算 | 理论距离 ($\Delta z$) |
| :--- | :--- | :--- |
| **L0 $\to$ L1** | $C_0=0.01$, $C_1=0.03$ | $0.02m$ |
| **L1 $\to$ L2** | $C_1=0.03$, $C_2=0.04+0.5(z_{snow}-0.04)$ | $0.5 z_{snow} - 0.01 \approx 0.5(z_{snow}-0.02) + 0.01$ |
| **L2 $\to$ Soil** | $C_2 \to \text{Interface} \to C_{soil}$ | Snow: $0.5(z_{snow}-0.04)$<br>Soil: $0.01m$ |

## 3. 代码修正方案

当前代码中的传导距离系统性偏大约 **2倍**。建议按下表修正：

### 3.1 变量修改对照表

| 物理过程 | 对应变量 | 当前代码 (Current) | **修正方案 (Proposed)** |
| :--- | :--- | :--- | :--- |
| **Top Implicit** | `solve_imp` ($z$) | `dz_snow_s12` (0.04) | **`dz_snow_s1` (0.02)** |
| **Flux L0$\to$L1** | `G_snow` (Dist) | `dz_snow_s12` (0.04) | **`dz_snow_s1` (0.02)** |
| **Flux L1$\to$L2** | `G_snow1` (Dist) | `z_snow - dz_snow_s1` | **`0.5 * (z_snow - dz_snow_s1) + 0.5 * dz_snow_s2`** |
| **Flux L2$\to$Soil** | `G_snow2` (Soil Dist) | `dz_soil_s0` (0.02) | **`0.5 * dz_soil_s0` (0.01)** |

### 3.2 推荐实现代码

```julia
# Case 3: 深雪 (>5cm) - 3层雪模型
dz_snow_s12 = dz_snow_s1 + dz_snow_s2

# 1. Top Layer (Implicit) - 使用 dz_snow_s1 作为传导距离
ΔM_snow = cp_ice * ρ_snow * dz_snow_s1 / Δt
T_snow = solve_imp(..., ra_g, dz_snow_s1, ...) # z changed to dz_snow_s1

# 2. Flux L0 -> L1
G_snow = κ_dry_snow * (T_snow - T_snow1_last) / dz_snow_s1

# 3. Flux L1 -> L2
dist_s1_s2 = 0.5 * (z_snow - dz_snow_s12) + 0.5 * dz_snow_s2 + 0.5 * dz_snow_s1 # (Simplified: 0.5*z_snow)
# 注: 严格推导距离为 0.5*z_snow - 0.01 + 0.5*0.02 = 0.5*z_snow
# 但考虑到节点位置, 推荐使用通用形式:
dist_s1_s2 = 0.5 * dz_snow_s2 + 0.5 * (z_snow - dz_snow_s12) # L1下半 + L2上半
G_snow1 = κ_dry_snow * (T_snow1_last - T_snow2_last) / dist_s1_s2

# 4. Flux L2 -> Soil
# 土壤侧热阻使用半层厚度 (0.5 * dz_soil_s0)
R_soil = (0.5 * dz_soil_s0) / κ_soil1
R_snow = (0.5 * (z_snow - dz_snow_s12)) / κ_dry_snow
G_snow2 = (T_snow2_last - T_soil_surf_last) / (R_snow + R_soil)
```
