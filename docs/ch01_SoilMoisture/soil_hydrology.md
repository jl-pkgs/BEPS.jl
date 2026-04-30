# 土壤水文过程参数说明

本文档介绍 BEPS.jl 中土壤水文模块的关键参数。

## 1. r_drainage (地表径流调节系数)

### 定义
`r_drainage` 是描述地表微地形对降水蓄留与排水能力的经验系数。它决定了在每一个计算步长中，未能下渗的“过剩水”中有多少比例被截留在地表积水洼地中。

### 物理逻辑
在 `src/Soil/UpdateSoilMoisture.jl` 中，地表积水深度 `z_water` 的更新逻辑如下：

```julia
# Ponded water after runoff.
soil.z_water = (soil.z_water / kstep + soil.r_rain_g - inf) * kstep * soil.r_drainage
```

对应的数学公式为：

$$ z_{water}^{t+1} = \left( \frac{z_{water}^{t}}{\Delta t} + P - I \right) \cdot \Delta t \cdot r_{drainage} $$

地表径流（Surface Runoff）则为剩余部分：

$$ R_{surface} = \left( \frac{z_{water}^{t}}{\Delta t} + P - I \right) \cdot \Delta t \cdot (1 - r_{drainage}) $$

其中：
- $z_{water}$: 地表积水深度 [m]
- $P$: 到达地表的降水速率 (`r_rain_g`) [m/s]
- $I$: 实际下渗速率 (`inf`) [m/s]，受限于最大入渗能力和积水量。
- $\Delta t$: 时间步长 (`kstep`) [s]
- $r_{drainage}$: 排水调节系数 [-]

### 物理意义与取值参考
- **物理意义**: 该参数反映了地形的“平坦程度”和“排水效率”。
- **取值范围**: 建议范围为 `[0.2, 0.7]`。
    - **低值 (如 0.2)**: 代表地形坡度大或排水系统通畅，积水很快转化为径流。
    - **高值 (如 0.7)**: 代表地形平坦、洼地多，水分容易滞留。
- **默认值**: 0.5。

### 对模型的影响
- **对含水量的影响**: 较大的 `r_drainage` 会增加地表积水的停留时间，从而在降雨停止后提供更长时间的下渗驱动力，通常会提高深层土壤的含水量。
- **对产流的影响**: 它是控制地表径流峰值的重要参数。较小的 `r_drainage` 会导致径流过程线更加陡峭（响应更快）。

## 2. 土壤水分传输方程

BEPS 采用基于 Richards 方程的分层模型来计算土壤水分的垂直运移。

### 2.1 达西定律离散化
第 $i$ 层与第 $i+1$ 层之间的水分通量 $Q_i$ [m/s] 计算公式如下：

$$ Q_i = K_{avg} \cdot \left( \frac{2(\psi_{i+1} - \psi_i)}{dz_i + dz_{i+1}} + 1 \right) $$

其中：
- $K_{avg}$: 层间平均导水率，代码中使用加权平均计算 (`KK`)。
- $\psi$: 土壤基质势 (Matric Potential) [m]。
- $dz$: 土壤层厚度 [m]。
- $+1$: 代表重力势梯度的贡献（向下为正）。

### 2.2 水分平衡更新
各层土壤含水量 $\theta_i$ 的随时间变化率由流入通量、流出通量和根系吸水共同决定：

$$ \frac{\partial \theta_i}{\partial t} = \frac{Q_{i-1} - Q_i - S_i}{dz_i} $$

离散化更新公式：

$$ \theta_i^{t+1} = \theta_i^t + \left( Q_{i-1} - Q_i - E_{tt,i} \right) \cdot \frac{\Delta t}{dz_i} $$

其中：
- $Q_{i-1}$: 上边界流入通量（对于第1层，$Q_0 = I$ 即入渗率）。
- $Q_i$: 下边界流出通量。
- $E_{tt,i}$: 该层的根系吸水速率 + 土壤蒸发速率 [m/s]。

## 3. r_root_decay (根系分布衰减率)

### 定义
描述根系随土壤深度增加而减少的指数衰减速率。

### 作用
用于计算各土壤层的 `f_root`（根系比例），直接影响：
1. **水分提取**: 决定了植被蒸腾从各层土壤吸水的比例。
2. **养分传输**: 影响根系相关的生化过程。

---
*更新日期: 2026-01-13*
