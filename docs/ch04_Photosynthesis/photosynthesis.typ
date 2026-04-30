#import "@local/modern-cug-report:0.1.3": *
#show: doc => template(doc, footer: "BEPS.jl", header: "")

#set par(first-line-indent: 2em)

#let pi = markhl.with(color: yellow)
#let pj = markhl.with(color: red)

#let LE = "LE"
#let RH = "RH"

= 1 植被光合

BEPS.jl 的光合模块基于 Farquhar-Ball-Berry 耦合框架，描述 CO#sub[2] 从大气扩散至叶绿体、
羧化反应产生光合产物、以及气孔动态响应的完整过程。


== 1.1 CO#sub[2] 扩散路径

CO#sub[2] 从大气（$c_a$）经边界层和气孔向胞间扩散，服从菲克扩散定律：

$
  c_s = c_a - A / g_b
$

$
  c_i = c_a - A (1/g_b + 1/g_s)
$

其中 $c_a$（ppm）为大气 CO#sub[2] 浓度，$c_s$ 为叶表面浓度，$c_i$ 为胞间浓度，
$g_b$（$mu$mol m#super[-2] s#super[-1]）为边界层导度，$g_s$ 为气孔导度，
$A$（$mu$mol m#super[-2] s#super[-1]）为净光合速率。

净光合与总光合的关系：

$
  A = A_g - R_d
$

$R_d$ 为暗呼吸速率。光照条件下暗呼吸按 Amthor 建议降低 40%：

$
  R_d = cases(
    0.4 dot R_(d,25) dot f_T (e_"rd") quad & "PPFD" > 5 mu"mol m"^(-2)"s"^(-1),
    R_(d,25) dot f_T (e_"rd") quad & "否则"
  )
$

其中 $R_(d,25) = 0.004657 dot V_(c"max",25)$。


== 1.2 Farquhar 光合模型

=== 1.2.1 统一形式

Rubisco 限制（$W_c$）和光能限制（$W_j$）均可写成同一结构：

$
  A_g = (a (c_i - Gamma)) / (e dot c_i + b)
$

#figure(
  table(
    columns: (auto, auto, auto, auto, auto),
    rows: (0.8cm, 0.9cm, 0.9cm),
    align: (horizon),
    table.header([*限制类型*], [$a$], [$e$], [$b$], [*物理含义*]),
    [Rubisco ($W_c$)], [$V_"cmax"$], [$1$], [$K_c (1 + O_2/K_o)$], [羧化酶底物饱和],
    [光能 ($W_j$)], [$J_x$], [$4$], [$8 Gamma$], [电子传递链限制],
  ),
  caption: [Farquhar 模型参数对照],
)

$Gamma$（$mu$mol mol#super[-1]）为 CO#sub[2] 补偿点（不含暗呼吸），由 Rubisco 特异性因子 $tau$ 决定：

$
  Gamma = (0.5 dot O_2) / tau times 1000
$

$K = K_c (1 + O_2 / K_o)$ 为表观 Michaelis-Menten 常数，综合了 CO#sub[2] 和 O#sub[2] 的竞争。

=== 1.2.2 极限行为分析

$
  c_i = Gamma => A_g = 0 quad ("CO"_2 "补偿点，羧化 = 光呼吸")
$
$
  c_i -> infinity => A_g -> a/e quad ("底物饱和，Rubisco 全部占满")
$

=== 1.2.3 光能限制项

电子传递速率 $J_x$ 由入射光（PPFD）和 $J_"max"$ 共同决定（Chen 1999, Eq. 6）：

$
  J_x = (J_"max" dot "PPFD") / ("PPFD" + 2.1 J_"max")
$

PPFD 由太阳短波辐射换算：$"PPFD" = 4.55 times 0.5 times R_"sn"$（$mu$mol m#super[-2] s#super[-1]）。

=== 1.2.4 蔗糖合成限制

Collatz 建议额外考虑蔗糖输出对光合的限制（第三限制）：

$
  A_"sucrose" = V_"cmax" / 2 - R_d
$

最终净光合取三者中最小值：

$
  A = min(W_c - R_d, quad W_j - R_d, quad A_"sucrose")
$


== 1.3 温度响应

=== 1.3.1 Arrhenius 函数（$K_c$, $K_o$, $tau$, $R_d$）

$
  f_T (e_"act") = exp((T_l - 25) dot e_"act" / (T_(K,25) dot R dot T_l))
$

其中 $T_l$（K）为叶温，$R = 8.314$ J mol#super[-1] K#super[-1] 为气体常数，
$e_"act"$ 为活化能（J mol#super[-1]）。

=== 1.3.2 Maxwell-Boltzmann 函数（$V_"cmax"$, $J_"max"$）

具有温度最适点的钟形温度响应（Harley & Tenhunen 1991）：

$
  f_T^"boltz" = (H_k dot exp(e_a (T_l - T_"opt") / (R T_"opt" T_l))) / (H_k - e_a (1 - exp(H_k (T_l - T_"opt") / (R T_"opt" T_l))))
$

其中 $H_k = 200000$ J mol#super[-1]，$T_"opt"$ 为最适温度。

$J_"max"$ 与 $V_"cmax"$ 的关系（Chen 1999, Eq. 7）：

$
  j_"mopt" = 2.39 dot V_(c"max",25) - 14.2
$


== 1.4 Ball-Berry 气孔导度模型

气孔对光合、湿度和 CO#sub[2] 的综合响应（Ball et al. 1987）：

$
  g_(s,w) = g_(0,w) + g_(1,w) dot "RH"_l dot beta_"soil" dot A / c_s
$

其中 $g_0$（最小导度）、$g_1$（斜率）为植被参数，$"RH"_l$ 为叶表面相对湿度，$beta_"soil"$ 为土壤水分胁迫因子（0–1）。$g_(s,c) = g_(s,w) / 1.6$（H#sub[2]O 与 CO#sub[2] 扩散比）。

叶表面相对湿度由叶温处饱和水汽压与潜热通量反推：

$
  rho_v = (LE / lambda) dot r_v + rho_a, quad
  e = rho_v T_l / 0.2165, quad
  "RH"_l = 1 - (e_s - e) / e_s
$


== 1.5 耦合方程的推导与求解

=== 1.5.1 二次方程推导（$g_s = g_0$，弱光或 $A < 0$）

$
  cases(
    // c_s = c_a - A 1 / g_b,
    g_s = g_0,
    c_i = c_a - A (1/g_b + 1/g_0),
    A_g = (a (c_i - Gamma)) / (e c_i + b),
    A = A_g - R_d
  )
$

当 $A < 0$ 时，气孔导度退化为纯截距 $g_s = g_0$。同时，将 $alpha = 1/g_b + 1/g_0$, $c_i = c_a - A alpha$ 代入 Farquhar 方程 $(A + R_d)(e c_i + b) = a(c_i - Gamma)$：

$
  (A + R_d)[e(c_a - A alpha) + b] = a(c_a - A alpha - Gamma)
$

左侧展开并移项，按 $A$ 的次幂整理：

$
  underbrace(-e alpha)_(alpha_2) A^2
  + underbrace((e c_a + b + alpha(a - e R_d)))_(beta_2) A
  + underbrace(R_d(e c_a + b) - a(c_a - Gamma))_(gamma_2) = 0
$

解为（取判别式 $Delta = beta_2^2 - 4 alpha_2 gamma_2 >= 0$ 时的根）：

$
  A = frac(-beta_2 + sqrt(Delta), 2 alpha_2)
$

注意 $alpha_2 = -e alpha < 0$，抛物线开口向下，上式给出较小的正根，与弱光条件吻合。

=== 1.5.2 三次方程推导（Ball-Berry 耦合）

方程组为（$c = g_1 dot RH dot beta_"soil"$，以 CO#sub[2] 导度计）：

$
  cases(
    c_s = c_a - A / g_b,
    g_s = g_0 + c A / c_s,
    c_i = c_s - A / g_s,
    A_g = (a (c_i - Gamma)) / (e c_i + b)\,,
    A = A_g - R_d
  )
$

策略：分别从 Farquhar 和 Ball-Berry 两条路径推出 $c_i$ 的表达式，令二者相等后消去 $c_i$，再代入 $c_s$，整理为关于 $A$ 的多项式。

*路径一：由 Farquhar 方程解出 $c_i$。*

$A_g = A + R_d = a(c_i - Gamma)/(e c_i + b)$，展开并解出 $c_i$：

$
  (A + R_d)(e c_i + b) = a c_i - a Gamma
  quad => quad
  c_i = frac(a Gamma + b(A + R_d), a - e(A + R_d))
$

*路径二：由气孔导度方程解出 $c_i$。*

将 $g_s = g_0 + c A / c_s$ 代入 $c_i = c_s - A/g_s$：

$
  c_i = c_s - frac(A c_s, g_0 c_s + c A) = frac(c_s (g_0 c_s + (c-1) A), g_0 c_s + c A)
$

*联立两路径，令 $c_i$ 相等：*

$
  frac(a Gamma + b(A + R_d), a - e(A + R_d))
  = frac(c_s [g_0 c_s + (c-1) A], g_0 c_s + c A)
$

交叉相乘，并代入 $c_s = c_a - A/g_b$：

$
  [pi(a Gamma + b(A + R_d))][g_0 c_s + c A] = [pj(a - e(A + R_d))] c_s [g_0 c_s + (c-1) A]
$ <eq_An3_main>

逐个整理各因子：

#box-red[
  *引入缩写。* 为避免后续展开过于冗长，记（下文都是常量）：
  $ M = a Gamma + b R_d, quad N = a - e R_d, quad theta = c g_b - g_0, quad alpha = 1 + g_0/g_b - c $
]

$ pi(a Gamma + b(A + R_d)) = (a Gamma + b R_d) + b A = M + b A $

$ pj(a - e(A + R_d)) = (a - e R_d) - e A = N - e A $

$
  g_0 c_s + c A = g_0(c_a - A/g_b) + c A = g_0 c_a + (c - g_0/g_b) A = bold(g_0 c_a + theta/g_b A)
$

$
  g_0 c_s + (c-1)A = g_0 c_a + ((c - 1) - g_0/g_b) A = g_0 c_a - alpha A
$

式#[@eq_An3_main]变为：

$
  (M + b A)(g_0 c_a + theta/g_b dot A) = (N - e A)(c_a - A/g_b)(g_0 c_a - alpha A)
$

*两侧同乘 $g_b$ 清除分母。*右侧的 $(c_a - A/g_b)(g_0 c_a - alpha A) dot g_b$ 展开为：

$
  & = (g_b c_a - A)(g_0 c_a - alpha A) \
  & = g_0 g_b c_a^2 - (alpha g_b + g_0) c_a A + alpha A^2
$

记 $gamma = g_0 g_b c_a^2$，$beta = -(alpha g_b + g_0)c_a = - c_a (g_b + 2g_0 - c g_b)$，则此二次式为 $gamma + beta A + alpha A^2$。同时 $g_0 c_a g_b = gamma/c_a$。

乘 $g_b$ 后的完整等式为：

$ (M + b A)(gamma/c_a + theta A) = (N - e A)(gamma + beta A + alpha A^2) $

*展开左右两端。*

左端：
$ (M + b A)(gamma/c_a + theta A) = M gamma/c_a + (M theta + b gamma/c_a) A + b theta A^2 $

右端：
$
  (N - e A)(gamma + beta A + alpha A^2)
  = N gamma + (N beta - e gamma) A + (N alpha - e beta) A^2 - e alpha A^3
$

*右减左令其为零，整理。*

$
  e alpha A^3 + (e beta + b theta - N alpha) A^2
  + (e gamma - N beta + M theta + b gamma/c_a) A
  + (M gamma/c_a - N gamma) = 0
$

*将 $M = a Gamma + b R_d$，$N = a - e R_d$ 代回，展开各系数。*

$A^2$ 系数：$e beta + b theta - N alpha = e beta + b theta - (a - e R_d) alpha = e beta + b theta - a alpha + e R_d alpha$

$A^1$ 系数：$e gamma - N beta + M theta + b gamma/c_a$
$= e gamma + b gamma/c_a - (a - e R_d) beta + (a Gamma + b R_d) theta$
$= e gamma + b gamma/c_a - a beta + a Gamma theta + e R_d beta + R_d b theta$

$A^0$ 系数：$M gamma/c_a - N gamma = (a Gamma + b R_d)gamma/c_a - (a - e R_d) gamma$
$= -a gamma + a Gamma gamma/c_a + e R_d gamma + R_d b gamma/c_a$

两侧除以 $e alpha$，得标准三次方程 $A^3 + p A^2 + q A + r = 0$：

$
  p &= frac(e beta + b theta - a alpha + e R_d alpha, e alpha) \
  q &= frac(e gamma + b gamma/c_a - a beta + a Gamma theta + e R_d beta + R_d b theta, e alpha) \
  r &= frac(-a gamma + a Gamma gamma/c_a + e R_d gamma + R_d b gamma/c_a, e alpha)
$

=== 1.5.3 三次方程的三角解法

令 $A = x - p/3$ 消去二次项，得压缩三次方程 $x^3 + alpha' x + beta' = 0$。
定义：

$
  Q = (p^2 - 3q)/9, quad U = (2p^3 - 9 p q + 27 r)/54
$

当 $U^2 <= Q^3$ 时方程有三个实根（Press, _Numerical Recipes_，三角解法）：

$
  psi = arccos(U / sqrt(Q^3))
$

$
  A_k = -2sqrt(Q) cos((psi + 2 k pi) / 3) - p/3, quad k = 0, 1, -1
$

*物理根的选取规则*（`findroot` 函数）：

#figure(
  table(
    columns: (auto, auto),
    align: (left, left),
    rows: (0.8cm, 0.8cm, 0.8cm, 0.8cm),
    table.header([*三根符号*], [*选取规则*]),
    [全为正根], [取最小正根（光合最弱的物理解）],
    [两负一正], [取唯一正根],
    [一负两正], [取中间正根],
  ),
  caption: [三次方程物理根选取规则],
)


== 1.6 计算流程

```
输入: T_leaf, Rsn, ea, gb_w, Vcmax25, β_soil, g0_w, g1_w, ci_old, LH_leaf, ca
  │
  ├─ 1. 温度响应: Kc, Ko, τ → Γ, K, Rd, Vcmax, Jmax
  ├─ 2. 光能: PPFD → Jx
  ├─ 3. 初始 ci 估算 Wc, Wj → 选择限制项 (a, b, e)
  ├─ 4. 叶面湿度: LE → RH_leaf → Ball-Berry 参数
  ├─ 5. 求解耦合方程:
  │     if Wj > Rd && Wc > Rd → solve_cubic (Ball-Berry)
  │     else                   → solve_quad  (gs = g0)
  │     An = min(An, j_sucrose); An = max(An, 0)
  └─ 6. 后处理: cs = ca - An/gb; gs_w; ci = cs - An/gs_c

输出: gs_w [s m⁻¹], An [μmol m⁻² s⁻¹], ci [ppm]
```

四类叶片（冠层/下木 × 阳叶/阴叶）分别调用 `photosynthesis_jl`，
阳叶/阴叶使用不同 $V_"cmax"$，冠层/下木使用不同边界层导度 $g_b$。
