#import "@local/modern-cug-report:0.1.3": *
#show: doc => template(doc, footer: "BEPS.jl", header: "")
#let MO = "MO"


= 1 Monin-Obukhov 相似理论（MOST）风速廓线推导

== 1.1 基本关系（公式1）

MOST 给出无量纲风速梯度与稳定度参数的关系（Bonan 2019, Eq 6.34）：

$ k (z - d) / u_star pdv(u, z) = phi_m ((z - d) / L_MO) $ <eq1>

其中 $k$ 为 von Kármán 常数，$d$ 为零平面位移高度，$u_star$ 为摩擦速度，$L_MO$ 为 Monin-Obukhov 长度。

== 1.2 推导过程

*第一步：整理 $dif u$*。由 @eq1 直接整理：

$ dif u = u_star / k dot (phi_m ((z-d) / L_MO)) / (z - d) dif z $

*第二步：无量纲换元*。令 $xi = (z - d) / L_MO$，则 $z - d = xi L_MO$，$dif z = L_MO dif xi$：

$
  dif u = u_star / k dot (phi_m (xi)) / (xi L_MO) dot L_MO dif xi
  = u_star / k dot (phi_m (xi)) / xi dif xi
$

*第三步：从动量粗糙度高度 $z_"0m"$ 到 $z$ 积分*。
如前所述，风速为零的绝对高度是 $z = d + z_"0m"$，边界条件为 $u(d + z_"0m") = 0$。

令 $xi = (z-d) / L_MO$，积分下限取 $xi_0 = z_"0m" / L_MO$，对应的绝对高度为：

$ xi_0 = z_"0m" / L_MO quad arrow.r.double quad z - d = z_"0m" quad arrow.r.double quad z = d + z_"0m" $

对 $dif u$ 两端从 $xi_0$ 到 $xi$ 积分：

$
  integral_(u(d + z_"0m"))^(u(z)) dif u
  = u_star / k integral_(xi_0)^xi (phi_m (x)) / x dif x
$

左端为 $u(z) - u(d + z_"0m") = u(z)$，故：

$ u(z) = u_star / k integral_(xi_0)^xi (phi_m (x)) / x dif x $

*第四步：拆分被积函数*。将$phi_m\/x$ 分解为 $1\/x$ 与 $(phi_m - 1)\/x$ 两部分：

$ u(z) = u_star / k [integral_(xi_0)^xi 1/x dif x + integral_(xi_0)^xi (phi_m (x) - 1) / x dif x] $

$ = u_star / k [ln xi/xi_0 - integral_(xi_0)^xi (1 - phi_m (x)) / x dif x] $

*第五步：引入稳定度修正函数 $Phi_m$*

$ Phi_m (xi) = integral_0^xi (1 - phi_m (x)) / x dif x $ <eq-Phi>

$ integral_(xi_0)^xi (1 - phi_m (x)) / x dif x = Phi_m (xi) - Phi_m (xi_0) $

代回，并注意 $xi\/xi_0 = (z-d)\/z_"0m"$：

== 1.3 积分结果（公式2）

$ u(z) = u_star / k [ln (z - d) / z_"0m" - Phi_m ((z - d) / L_MO) + Phi_m (z_"0m" / L_MO)] $ <eq2>

其中稳定度修正函数 $Phi_m$ 的定义见式@eq-Phi。

#block(
  fill: luma(230),
  inset: 8pt,
  radius: 4pt,
)[
  *说明：*$Phi_m (z_"0m" \/ L_MO)$ 项在大多数情况下量值很小，常省略；但保留该项可使推导在任意稳定度条件下严格成立。中性条件下 $phi_m equiv 1$，$Phi_m equiv 0$，式@eq2 退化为对数风速廓线。
]

= 2 空气动力学阻力 $r_a$ 的推导

== 2.1 MOST 温度廓线

类比动量传输，MOST 给出*热量*的通量-梯度关系：

$ (k (z - d)) / theta_star pdv(theta, z) = phi_h ((z - d) / L_MO) $ <eq-heat>

其中 $theta_star$ 为温度尺度，与感热通量 $H$ 的关系为：

$ H = -rho c_p u_star theta_star $ <eq-theta-star>

== 2.2 积分得温度廓线

对 @eq-heat 整理并积分（方法与风速廓线完全类比），从标量粗糙度高度 $z_"0h"$ 到观测高度 $z$：

$
  theta(z) - theta_s
  = theta_star / k
  [ln (z - d) / z_"0h" - Psi_h ((z-d) / L_MO) + Psi_h (z_"0h" / L_MO)]
$ <eq-T-profile>

其中 $theta_s = theta(d + z_"0h")$ 为地表温度，$Psi_h$ 为热量稳定度修正函数（定义形式同 @eq-Phi）。

== 2.3 定义空气动力学阻力

感热通量的体积传输公式定义 $r_a$（单位 s/m）：

$ H = rho c_p (theta_s - theta(z)) / r_a $ <eq-ra-def>

由 @eq-theta-star 得 $theta_star = -H \/ (rho c_p u_star)$，代入 @eq-T-profile：

$
  theta_s - theta(z)
  = - theta_star / k [ln (z-d) / z_"0h" - Psi_h]
  = H / (rho c_p u_star k) [ln (z-d) / z_"0h" - Psi_h]
$

与 @eq-ra-def 比较，得：

$ r_a = 1 / (k u_star) [ln (z - d) / z_"0h" - Psi_h ((z-d) / L_MO) + Psi_h (z_"0h" / L_MO)] $ <eq-ra>

通常省略量值很小的 $Psi_h (z_"0h" \/ L_MO)$ 项，简化为：

$ r_a approx 1 / (k u_star) [ln (z - d) / z_"0h" - Psi_h ((z-d) / L_MO)] $

#block(
  fill: luma(230),
  inset: 8pt,
  radius: 4pt,
)[
  *注意：* 公式中的粗糙度长度是*热量粗糙度长度* $z_"0h"$，而非*动量粗糙度长度* $z_"0m"$，两者通过 $k B^(-1)$ 参数关联：$ln(z_"0m" \/ z_"0h") = k B^(-1)$，典型值 $z_"0h" approx 0.1 z_"0m"$。若误用 $z_"0m"$ 代替 $z_"0h"$，会高估 $r_a$（低估感热通量）。
]

= 3 摩擦速度 $u_star$ 的计算

== 3.1 由动量廓线反解

$u_star$ 由@eq2（动量廓线）直接反解。@eq2 的完整形式为：

$ u(z) = u_star / k [ln (z-d) / z_"0m" - Psi_m ((z-d) / L_MO) + Psi_m (z_"0m" / L_MO)] $

反解 $u_star$：

$ u_star = (k dot u(z)) / (ln (z-d) / z_"0m" - Psi_m ((z-d) / L_MO) + Psi_m (z_"0m" / L_MO)) $

由于 $z_"0m" << z - d$，故 $z_"0m" / L_MO << (z-d) / L_MO$，稳定度修正函数在该点的值接近于零：

$ Psi_m (z_"0m" / L_MO) approx 0 $

因此简化为：

$ u_star approx (k dot u(z)) / (ln (z-d) / z_"0m" - Psi_m ((z-d) / L_MO)) $ <eq-ustar>

其中积分下边界为*动量粗糙度长度* $z_"0m"$（而非热量粗糙度长度 $z_"0h"$），$Psi_m$ 为动量稳定度修正函数。

== 3.2 循环依赖与中性近似

@eq-ustar 存在循环依赖：计算 $u_star$ 需要 $L_MO$，而 $L_MO$ 又依赖 $u_star$：

$ L_MO = -(rho c_p (T_"air" + 273.15) u_star^3) / (k g H) $

实践中通常采用以下方案之一：

+ *中性近似（一次计算）：* 令 $Psi_m = 0$，得初步估算值：
  $ u_star approx (k dot u(z)) / ln((z-d) \/ z_"0m") $
  随后用此 $u_star$ 计算 $L_MO$ 和 $xi$，再修正 $r_a$。

+ *迭代法（更精确）：* 以中性近似值为初值，迭代 2–3 次直至收敛。

== 3.3 计算流程小结

$ u_star^((0)) = (k u) / ln((z-d) \/ z_"0m") arrow.r L_MO arrow.r xi = (z-d)/L_MO arrow.r u_star^((1)) = (k u) / (ln((z-d) \/ z_"0m") - Psi_m (xi)) $

最终 $r_a$ 使用稳定后的 $u_star$ 和热量粗糙度 $z_"0h"$ 计算（@eq-ra）。

#block(
  fill: luma(230),
  inset: 8pt,
  radius: 4pt,
)[
  *重要区别：*
  - $u_star$ 计算用 $z_"0m"$（动量廓线积分下边界）
  - $r_a$ 计算用 $z_"0h"$（热量/标量廓线积分下边界）
  - 典型关系：$z_"0h" approx 0.1 z_"0m"$（$k B^(-1) approx 2.3$）
  
  混用两者是常见错误，会同时导致 $u_star$ 偏高和 $r_a$ 偏低。
]

= 4 冠层内部空气动力学阻力 $r_a$

== 4.1 物理假设与控制方程

基于一阶闭合湍流扩散理论（K-theory），冠层内的湍流通量 $F$ 与标量浓度 $C$ 的垂直梯度成正比：
$ F = -K_H(z) (partial C) / (partial z) $

由于冠层内部枝叶的机械阻挡，湍流被耗散。假定湍流热扩散系数 $K_H(z)$ 与风速 $u(z)$ 随深度呈指数衰减：
$
  K_H(z) & = K_h exp(-gamma (1 - z/h)) \
    u(z) & = u_h exp(-gamma (1 - z/h))
$

其中，$h$ 为冠层顶部高度；$K_h$ 为冠层顶部的湍流扩散系数；$gamma$ 为湍流消光系数（Extinction coefficient），反映枝叶的阻挡强度。

== 4.2 空气动力学阻力的通用积分

两点 $z_1$ 和 $z_2$（且 $z_2 > z_1$）之间的空气动力学阻力 $r_a$ 定义为湍流扩散系数倒数在垂直方向上的积分：
$ r_a(z_1, z_2) = integral_(z_1)^(z_2) 1 / (K_H(z)) dif z $

将指数衰减方程代入：
$ r_a(z_1, z_2) = 1 / K_h integral_(z_1)^(z_2) exp(gamma (1 - z/h)) dif z $

引入换元法，令 $x = 1 - z/h$，则 $dif z = -h dif x$。积分上下限相应变为 $x_1 = 1 - z_1/h$ 和 $x_2 = 1 - z_2/h$：
$ r_a = h / K_h integral_(x_2)^(x_1) exp(gamma x) dif x $

完成指数函数的积分计算，并替换回 $z$ 的表达式，得到阻力的通用计算公式：
$
  r_a(z_1, z_2) & = h / (gamma K_h) [exp(gamma x)]_(x_2)^(x_1) \
                & = h / (gamma K_h) [exp(gamma (1 - z_1/h)) - exp(gamma (1 - z_2/h))]
$

== 4.3 分层阻力的精确计算（$gamma$ 相同场景）

=== 4.3.1 林下植被到冠层顶部的阻力 ($r_(a,u)$)
积分区间为 $z_1 = h_u$ 到 $z_2 = h$。代入通用公式：
$
  r_(a,u) & = h / (gamma K_h) [exp(gamma (1 - h_u/h)) - exp(gamma (1 - h/h))] \
          & = h / (gamma K_h) [exp(gamma (1 - h_u/h)) - 1]
$

=== 4.3.2 地面到林下植被的阻力 ($r_(a,g)$)
积分区间为地面 $z_1 = 0$ 到 $z_2 = h_u$。代入通用公式：
$
  r_(a,g) & = h / (gamma K_h) [exp(gamma (1 - 0/h)) - exp(gamma (1 - h_u/h))] \
          & = h / (gamma K_h) [exp(gamma) - exp(gamma (1 - h_u/h))]
$

== 4.4 BEPS 模型中的实现细节（分层 $gamma$ 场景）

在真实的复杂冠层中，流体动量（风速）的机械消耗与标量热量（湍流）的扩散衰减机制不同。为保证严谨性，模型中对两者的消光系数（$gamma$）进行了物理层面的解耦：
- $gamma_(o,m)$：主林冠层（Overstory）的*动量*消光系数，用于计算风速的衰减。
- $gamma_(o,h)$：主林冠层（Overstory）的*热量/标量*消光系数，用于计算主林冠下部及树干对湍流扩散系数的衰减与空气动力学阻力。
- $gamma_(u,m)$：林下植被层（Understory）的*动量*消光系数，用于计算林下空间风速的衰减阻滞。
- $gamma_(g,h)$：极近地表层（Ground）的*热量/标量*衰减常数，用于计算林下层底部到地表的极近地表层空气阻力。

=== 4.4.1 有效风速与叶片边界层阻力 ($r_b$)

风速在冠层内的衰减严格遵循动量消光系数。主林冠内的有效风速 $u_o$ 及林下有效风速 $u_u$ 的衰减路径为：
$
  u_o & = u_h exp(-gamma_(o,m) (1 - d/h)) \
  u_u & = u_"hu" exp(-gamma_(u,m) (1 - 0.8 h_u/h_u)), quad u_"hu" = u_h exp(-gamma_(o,m) (1 - h_u/h))
$

#box-blue[
  *有效风速的物理意义* \
  代码中计算衰减到零平面位移高度 $d$ 处的风速 $u_o$，是为了提取出一个*“最能代表广大叶片实际感受到的平均微风速度”*。用这个有效风速去算出来的努塞尔数（Nusselt number）和叶片边界层阻力 $r_b$，才最符合整棵树的宏观统计规律。
]

=== 4.4.2 林冠层空气动力学阻力 ($r_(a,u)$)
积分区间从林下植被顶部 $h_u$ 到大树顶部 $h$。此空间内的湍流受主林冠的热量消光系数 $gamma_(o,h)$ 控制：
$ r_(a,u) = h / (gamma_(o,h) K_(h,o)) [exp(gamma_(o,h) (1 - h_u/h)) - 1] $

=== 4.4.3 底层空间空气动力学阻力 ($r_(a,g)$)

+ *湍流强度的连续性交接* \
  在交界面（$z = h_u$）处，上下层湍流扩散系数必须保持物理连续。由于流体刚刚穿过主林冠空间，底层积分的基准 $K_(h,u)$ 必须由上层的*热量消光方程*向下推导求得：
  $ K_(h,u) = K_(h,o) exp(-gamma_(o,h) (1 - h_u/h)) $

+ *底层阻力积分* \
  在林下空间（$0$ 到 $h_u$），衰减基准切换为 $h_u$ 和 $K_(h,u)$，受底层消光系数 $gamma_(g,h)$ 控制。此时的剖面方程为 $K_H (z) = K_(h,u) exp(-gamma_(g,h) (1 - z/h_u))$。对此方程积分：
$
  r_(a,g) & = integral_0^(h_u) 1 / (K_(h,u) exp(-gamma_(g,h) (1 - z/h_u))) dif z \
          & = h_u / (gamma_(g,h) K_(h,u)) [exp(gamma_(g,h) (1 - 0/h_u)) - exp(gamma_(g,h) (1 - h_u/h_u))] \
          & = h_u / (gamma_(g,h) K_(h,u)) [exp(gamma_(g,h)) - 1]
$

#box-red[
  *物理拓扑提醒：阻力网络的并联属性* \
  地表通量与各层叶片通量通过各自的阻力（$r_(a,g)$, $r_(b,u)$, $r_(b,o)$）*并联*注入冠层内部的空气节点，最终统一经由 $r_(a,o)$ 排入大气层。各分层阻力代表独立的传输通道，*严禁将各层空气阻力进行简单的电路串联累加*（即 $r_(a,g) eq.not r_(a,g) + r_(a,u) + r_(a,o)$），否则会破坏能量守恒定律并导致底层通量极度失真。
]

// #pagebreak()

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
