#import "@preview/modern-cug-report:0.1.3": *
#show: doc => template(doc, footer: "BEPS.jl", header: "")
#let MO = "MO"

#v(2em)
= 5 V1 与 V2 实现对比

本模型经历了两个主要实现版本（`aerodynamic_conductance.jl` 与 `aerodynamic_conductance_V2.jl`），V2 在物理一致性上有系统性提升。

== 5.1 主要改进点汇总

#figure(
  table(
    columns: (auto, 1fr, 1fr),
    align: (left, left, left),
    table.header([*项目*], [*V1（旧版）*], [*V2（新版）*]),
    [Monin-Obukhov 长度 $L$],
    [$L = -(k g H) / (rho c_p T u_star^3)$（公式倒置）],
    [$L = -(rho c_p T u_star^3) / (k g H)$（正确 MOST 公式）],
    
    [热量粗糙度 $z_"0h"$],
    [与动量粗糙度混用：$z_0 = 0.08h$],
    [显式区分：$z_"0m" = 0.08h$，$z_"0h" = 0.1 z_"0m"$（$k B^(-1) approx 2.3$）],
    
    [$r_(a,o)$ 稳定性修正],
    [$1/(k u_star) (ln((z-d)\/z_0) + n (z-d) L)$（加法形式，错误）],
    [$1/(k u_star) (ln((z-d)\/z_"0h") - Psi_h)$（正确 MOST 积分形式）],
    
    [消光系数 $gamma$],
    [所有层共用同一 $gamma$],
    [动量/热量解耦：$gamma_(o,m)$、$gamma_(o,h)$、$gamma_(u,m)$、$gamma_(g,h)$],
    
    [$K_H$ 连续性],
    [林下阻力直接使用 $K_(h,o)$（冠层顶部值）],
    [$K_(h,u) = K_(h,o) exp(-gamma_(o,h)(1-h_u\/h))$（物理连续交接）],
    
    [$r_(a,g)$ 公式],
    [$h / (gamma K_(h,o)) [exp(gamma) - exp(gamma (1-h_u\/h))]$],
    [$h_u / (gamma_(g,h) K_(h,u)) [exp(gamma_(g,h)) - 1]$（以 $h_u$ 为基准，更简洁）],
  ),
  caption: [V1 与 V2 关键差异对比],
)

== 5.2 关键公式对比

=== Monin-Obukhov 长度（最重要的改正）

V1 与 V2 中 L 的计算互为倒数关系（$L_"V1" dot L_"V2" = 1$）：

$
  L_"V1" = -(k g H) / (rho c_p T u_star^3) quad [upright("量纲: m"^(-1))]
  quad arrow.r quad
  L_"V2" = -(rho c_p T u_star^3) / (k g H) quad [upright("量纲: m")]
$

$L_"V1"$ 量纲为 m⁻¹，与 Monin-Obukhov 长度的物理单位（m）不符，属于量纲错误。两者之间满足 $L_"V1" = 1 \/ L_"V2"$。

然而由于 V1 在计算稳定度参数时使用乘法 $xi_"V1" = (z-d) dot L_"V1"$，而 V2 使用除法 $xi_"V2" = (z-d) \/ L_"V2"$，两者恰好相等：

$
  xi_"V1" = (z-d) dot L_"V1" = (z-d) dot frac(1, L_"V2") = frac(z-d, L_"V2") = xi_"V2"
$

因此 *$xi$ 的数值计算结果一致*，V1 的 L 公式错误在稳定度参数层面被"抵消"。但在 V1 的 $r_(a,o)$ 公式中使用了加法形式的稳定性修正 $+n dot xi$，而非正确的积分形式 $-Psi_h$：

$
  r_(a,o)^"V1" = 1 / (k u_star) [ln((z-d) / z_"0m") + n xi]
  quad arrow.r quad
  r_(a,o)^"V2" = 1 / (k u_star) [ln((z-d) / z_"0h") - Psi_h (xi)]
$

=== 热量粗糙度长度

$r_(a,o)$ 的推导（见 @eq-ra）积分下边界为热量粗糙度长度 $z_"0h"$，而非动量粗糙度长度 $z_"0m"$。V1 误用 $z_"0m"$，而 V2 正确引入：

$
  z_"0h" = 0.1 z_"0m" = 0.008 h
  quad arrow.r quad
  ln((z-d) / z_"0h") = ln((z-d) / z_"0m") + ln(z_"0m" / z_"0h") = ln((z-d) / z_"0m") + k B^(-1)
$

其中 $k B^(-1) approx 2.3$。这使得 V2 的 $r_(a,o)$ 系统性地高于 V1 约 $k B^(-1) / (k u_star) approx 5.75 / u_star$ s/m。

== 5.3 V2 对感热通量模拟精度的影响

=== 5.3.1 感热通量传播链

感热通量从气动阻力到最终输出，经历以下传播链：

$
  r_(a,o) arrow.r G_(a,o) = 1\/r_(a,o) arrow.r G_h = 1\/(1\/G_(a,o) + 0.5\/G_(b,o)) arrow.r T_c arrow.r "SH" = rho_a c_p G_h (T_c - T_a)
$

其中叶温方程（忽略气孔水汽项的简化形式）为：

$
  T_c = T_a + frac(R_n - "VPD" dot rho_a c_p G_w \/ gamma, rho_a c_p (G_h + Delta G_w \/ gamma))
$ <eq-Tc>

V2 相对 V1 使 $r_(a,o)$ 增大（典型中性条件增大约 15 s/m，见 @tab-v1v2-compare），导致 $G_h$ 减小。

#figure(
  table(
    columns: (auto, 1fr, 1fr, 1fr),
    align: (left, center, center, center),
    table.header([*变量*], [*V1*], [*V2*], [*差异*]),
    [$r_(a,o)$ [s/m]], [6.87], [22.28], [+15.41],
    [$G_(a,o)$ [m/s]], [0.1455], [0.0449], [-69%],
    [$G_h$ [m/s]], [0.0276], [0.0150], [-46%],
    [$T_c$ [°C]], [22.3], [24.6], [+2.3 K],
    [SH [W/m²]], [64], [70], [+10%],
  ),
  caption: [典型条件（$h=20$ m, $u=3$ m/s, $H_"init"=100$ W/m², $R_n=300$ W/m²）下 V1 与 V2 的差异，],
) <tab-v1v2-compare>

=== 5.3.2 Bowen 比效应

@eq-Tc 显示，SH 对 $G_h$ 的敏感性取决于 $G_h$ 与 $Delta G_w\/gamma$ 的相对大小（即 Bowen 比 $beta$）。定义：

$
  beta = "SH" \/ "LE" = G_h \/ (Delta G_w \/ gamma)
$

- *干旱/低蒸腾条件*（$G_w → 0$）：$T_c ≈ T_a + R_n \/ (rho_a c_p G_h)$，SH 与 $G_h$ 近线性负相关。V2 的 $G_h$ 减小导致 SH 非单调变化，改进最显著。
- *湿润/高蒸腾条件*（$G_w$ 大）：$T_c - T_a$ 小，SH 对 $G_h$ 不敏感，V1/V2 差异减弱。

=== 5.3.3 逆向诊断方案

==== 方法 A：直接气动阻力反演（需冠层温度 LST）

若有红外测温 $T_c$ 和涡度协方差 SH 观测，可反演叶片总传热阻力 $r_H = 1\/G_h$：

$
  r_H = frac(rho_a c_p (T_c - T_a), "SH"_"obs") = 1\/G_h = r_(a,o) + 0.5 r_(b,o)
$ <eq-ra-inv>

进而分离气动阻力：$r_(a,o)^"EC" = r_H - 0.5 r_(b,o)$（$r_{b,o}$ 由叶片边界层模型提供）。将 $r_(a,o)^"EC"$ 与 V1/V2 预测值对比，可直接评估两版本精度（见代码 `ra_from_flux`）。

==== 方法 B：Penman-Monteith 逆推（仅需 LE 观测）

不需要 LST，从能量平衡逆推冠层导度：

$
  G_c = frac(gamma dot "LE"_"obs" dot G_a, Delta(R_n - G) + rho_a c_p "VPD" dot G_a - (Delta + gamma) "LE"_"obs")
$ <eq-Gc-inv>

$G_a = 1\/r_(a,o)$ 分别取 V1/V2 值，对比逆推得到的 $G_c$ 与模型 $G_w$ 之差，判断哪个版本更接近生理合理值（见代码 `Gc_penman`）。

==== 分稳定度区间统计

利用 `stability_class(ξ)` 函数将逐小时数据分为三类（$xi<-0.5$：不稳定；$-0.5~0.5$：中性；$>0.5$：稳定），分别计算 V1/V2 的 SH RMSE，量化改进效果。完整示例见 `examples/example_diagnose_ra.jl`。

#v(2em)
= 附录 // <!-- omit in toc -->
// =

== A1. 两个高度参数的物理含义 // <!-- omit in toc -->

*零平面位移高度 $d$（zero-plane displacement height）*

对于裸地，风速在地面（$z=0$）处降为零；但对于有植被覆盖的下垫面，气流受冠层阻挡，风速廓线的"有效零点"被抬升至冠层内某处。$d$ 就是这个等效零点距地面的高度。引入 $d$ 后，公式中所有高度均以 $z - d$ 表示，即相对于等效零点的距离，而非相对于地面。

*动量粗糙度长度 $z_"0m"$（momentum roughness length）*

将对数风速廓线向下外推，风速降为零时对应的高度即为 $z_"0m"$。这并非真实的物理高度，而是衡量地表对动量传输阻力的参数——$z_"0m"$ 越大，地表越粗糙，对气流阻力越强。以绝对高度 $z = d + z_"0m"$ 处作为积分下边界，即 $u(d + z_"0m") = 0$。

#figure(
  table(
    columns: (auto, 10cm, 1fr),
    align: (horizon, horizon, horizon),
    table.header([*参数*], [*物理意义*], [*典型估算*]),
    [$d$], [风速廓线等效零点相对地面的高度], [$approx 2\/3 h$],
    [$z_"0m"$], [地表粗糙度的等效长度，风速外推为零处], [$approx 0.1 h$],
  ),
  caption: [$d$ 与 $z_"0m"$ 的含义对比（$h$ 为植被冠层高度）],
)

*两者的几何关系：* $z_"0m"$ 是在以 $d$ 为原点的坐标系（即 $z-d$ 坐标）中量取的，因此风速为零的绝对高度为：

$ z_"ref" = d + z_"0m" $

即 $d$ 确定坐标原点，$z_"0m"$ 在此基础上给出积分下边界的位置。
