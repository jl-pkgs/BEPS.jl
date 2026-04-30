# 1 T_weighted 公式推导与代码审查

本文档基于能量守恒定律与傅里叶热传导定律，推导 `surface_temperature.jl` 中雪-土混合层温度 `T_weighted` 的理论公式，并检查代码实现的正确性。

## 1.1 物理模型与符号定义

为了统一描述，我们建立如下垂直一维坐标系与物理量定义：

### 1.1.1 坐标系与符号约定

- **坐标轴 (**$z$**)**：定义**垂直向上为正方向**。

  - 界面位置：$z = 0$。
  - 雪层节点：$z = z_{snow}$ ($z > 0$)。
  - 土壤节点：$z = -z_{soil}$ ($z < 0$，其中 $z_{soil}$ 为正值的深度距离)。

- **热通量 (**$G$**)**：定义**沿** $z$ **轴正方向（向上）为正**。

  - 遵循统一的傅里叶定律矢量形式：

    $$G = - \kappa \frac{\partial T}{\partial z}$$

  - 注：若 $\frac{\partial T}{\partial z} < 0$（温度随高度降低），则 $G > 0$（热量向上输送）。这符合“向下通量为负”的定义。

### 1.1.2 变量定义

| 物理量                    | 代码变量               | 含义                             |
| ------------------------- | ---------------------- | -------------------------------- |
| $T_w$                     | `T_weighted`           | 待求的界面混合层温度 ($z=0$)     |
| $T_{old}$                 | `T_soil_surf_last`     | 上一时刻界面温度                 |
| $T_{snow}$                | `T_snow`               | 上方雪层温度 ($z > 0$)           |
| $T_{soil}$                | `T_soil1_last`         | 下方土壤温度 ($z < 0$)           |
| $\kappa_{snow}, z_{snow}$ | `κ_dry_snow`, `z_snow` | 雪层热导率与层厚                 |
| $\kappa_{soil}, z_{soil}$ | `κ_soil1`, `z_soil1`   | 土壤热导率与层厚                 |
| $M_{surf}$                | `ΔM_soil1`             | 界面层热惯性系数 (热容/时间步长) |

## 1.2 能量守恒方程

对位于 $z=0$ 处的界面薄层控制体应用能量守恒定律。内能的变化率等于热通量的负散度（即净通量辐合）。

$$
M_{surf} (T_w - T_{old}) = -\nabla \cdot G \approx - G_{in} + G_{out}
$$

由于向下为负，$-G_{in}$为进来能量的绝对值；公式右侧表达的是，进来的能量减去出去的能量。

在垂直方向上，这等价于**从下方流入的通量**减去**从上方流出的通量**：
$$
M_{surf} (T_w - T_{old}) = - G|_{z=0^+} + G|_{z=0^-} 
$$

其中：

- $G|_{z=0^+}$ 为界面上侧通量（界面 $\to$ 雪层）。
- $G|_{z=0^-}$ 为界面下侧通量（土壤 $\to$ 界面）。

## 1.3 公式推导

### 1.3.1 离散化热通量

直接应用傅里叶定律 $G = - \kappa \frac{\partial T}{\partial z}$：

**1. 界面上侧 (**$z=0^+$**)**

$$
G|_{z=0^+} = - \kappa_{snow} \frac{T_{snow} - T_w}{z_{snow} - 0} = \frac{\kappa_{snow}}{z_{snow}} (T_w - T_{snow})
$$

**2. 界面下侧 (**$z=0^-$**)**

$$
G|_{z=0^-} = - \kappa_{soil} \frac{T_w - T_{soil}}{0 - (-z_{soil})} = \frac{\kappa_{soil}}{z_{soil}} (T_{soil} - T_w)
$$

### 1.3.2 代入守恒方程求解

将通量表达式代入能量守恒方程：

$$
M_{surf} (T_w - T_{old}) = - \left[ \frac{\kappa_{snow}}{z_{snow}} (T_w - T_{snow}) \right] + \left[ \frac{\kappa_{soil}}{z_{soil}} (T_{soil} - T_w) \right] 
$$

展开各项：

$$
M_{surf} T_w - M_{surf} T_{old} = \frac{\kappa_{soil}}{z_{soil}} T_{soil} - \frac{\kappa_{soil}}{z_{soil}} T_w - \frac{\kappa_{snow}}{z_{snow}} T_w + \frac{\kappa_{snow}}{z_{snow}} T_{snow}
$$

将含有未知量 $T_w$ 的项移至左侧，已知量移至右侧：

$$
T_w \left( M_{surf} + \frac{\kappa_{soil}}{z_{soil}} + \frac{\kappa_{snow}}{z_{snow}} \right) = \frac{\kappa_{soil}}{z_{soil}} T_{soil} + \frac{\kappa_{snow}}{z_{snow}} T_{snow} + M_{surf} T_{old}
$$

最终得到 $T_w$ 的解析解：

$$
T_{weighted} = \frac{ \frac{\kappa_{soil}}{z_{soil}} T_{soil} + \frac{\kappa_{snow}}{z_{snow}} T_{snow} + M_{surf} T_{old} }{ \frac{\kappa_{soil}}{z_{soil}} + \frac{\kappa_{snow}}{z_{snow}} + M_{surf} }
$$

## 1.4 代码审查与结论

### 1.4.1 代码实现对比

检查 `surface_temperature.jl` 中 Case 2 (Line 60-63) 的实现：

```
# 代码中的分子实现
Numerator = (κ_soil1 * T_soil1_last / z_soil1) + 
            (T_snow * κ_dry_snow) +           # <--- 关注此项
            (ΔM_soil1 * T_soil_surf_last)
```

对应的数学形式为：

$$\text{Numerator}_{code} = \frac{\kappa_{soil}}{z_{soil}} T_{soil} + \kappa_{snow} T_{snow} + M_{surf} T_{old}$$

==第二项漏了分子$z_{snow}$==

### 1.4.2 错误分析

对比推导结果与代码实现，可以发现**雪层热通量项存在量纲错误**。

- **理论项**：$\displaystyle \frac{\kappa_{snow}}{z_{snow}} T_{snow}$
  - 物理量纲：热导度 $\times$ 温度 $[W \cdot m^{-2} \cdot K^{-1}] \cdot [K] = [W \cdot m^{-2}]$ (热通量)
- **代码项**：$\displaystyle \kappa_{snow} T_{snow}$
  - 物理量纲：热导率 $\times$ 温度 $[W \cdot m^{-1} \cdot K^{-1}] \cdot [K] = [W \cdot m^{-1}]$ (**错误**)
