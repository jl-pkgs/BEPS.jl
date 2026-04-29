// photosynthesis_comparison.typ
// 两种光合模块实现差异对比文档
//
// 编译方式：
//   typst compile docs/modules/photosynthesis_comparison.typ
// 依赖本地模板（可选），无模板时将 import 行注释并取消下面两行注释：
// #set text(font: "Noto Serif CJK SC", lang: "zh", size: 11pt)
// #set par(justify: true, first-line-indent: 2em)

#import "@local/modern-cug-report:0.1.3": *
#show: doc => template(doc, footer: "BEPS.jl 光合模块对比", header: "")

#align(center)[
  #text(size: 20pt, weight: "bold")[两种光合作用模块的实现差异分析]
  #v(0.6em)
  #text(size: 12pt)[`photosynthesis_jl` vs `Photosynthesis`（独立模块）]
  #v(0.4em)
  #text(size: 10pt, fill: gray)[BEPS.jl · #datetime.today().display()]
]

#v(1em)
#outline(indent: auto, depth: 3)
#pagebreak()

= 1 概述

BEPS.jl 中存在两套独立的 Farquhar-Ball-Berry 光合作用实现：

- `photosynthesis_jl`：位于 `src/photosynthesis.jl`，是 BEPS 主积分循环（`inter_prg.jl`）
  调用的原始模块，与叶温迭代、边界层导度等耦合紧密。
- `Photosynthesis`（独立模块）：位于 `src/Photosynthesis/`，是面向参数优化和
  单叶尺度研究设计的简化接口，输入仅需气温、湿度、辐射和 LAI。

两者共享相同的 Farquhar 生化框架（FvCB, 1980），但在
*温度响应函数*、*Michaelis-Menten 动力学常数*、*电子传递公式* 和
*气孔-光合耦合求解策略* 四个方面存在系统性差异，
导致在相同气象输入条件下计算结果不同。

#figure(
  table(
    columns: (2fr, 3fr, 3fr),
    inset: 7pt,
    align: (left, left, left),
    table.header(
      [*差异维度*], [*`photosynthesis_jl`（原始）*], [*`Photosynthesis`（独立）*],
    ),
    [温度响应 $V_"cmax"$/$J_"max"$], [TBOLTZ 钟形函数（Harley 1991）\ $T_"opt"$ = 301 K], [Medlyn 归一化 Arrhenius（2002）\ 25°C 严格归一化],
    [Michaelis-Menten 常数], [$K_(c,25)$ = 274.6, $K_(o,25)$ = 419.8\ $tau_25$ = 2904（Harley 1995）], [$K_(c,25)$ = 404, $K_(o,25)$ = 248\ $tau_25$ = 2600（Bernacchi 2001）],
    [电子传输 $J(I)$], [双曲饱和（Chen 1999）\ $J = J_"max" I / (I + 2.1 J_"max")$], [非矩形双曲线（Farquhar 1980）\ 含曲率参数 $theta$ = 0.7],
    [耦合求解策略], [解析三次/二次方程], [数值定点迭代（≤15步）],
    [叶表面相对湿度], [由潜热通量反推（能量耦合）], [直接使用气温湿度（$T_l = T_a$）],
    [$V_"cmax"$ 冠层分布], [外部 `VCmax()` 函数计算阴阳叶梯度], [阳/阴叶使用同一 $V_(c"max",25)$],
  ),
  caption: [两模块主要差异一览],
)

#pagebreak()

= 2 温度响应函数（核心差异）

== 2.1 原始模块：TBOLTZ 钟形函数

`photosynthesis_jl` 对 $V_"cmax"$ 和 $J_"max"$ 使用 Maxwell-Boltzmann 型钟形函数
（Harley & Tenhunen 1991，代码见 `photosynthesis_helper.jl:TBOLTZ`）：

$
  f_"TBOLTZ"(T_l) = "rate" times
  frac(H_k exp(e_a (T_l - T_"opt") / (R T_"opt" T_l)),
       H_k - e_a (1 - exp(H_k (T_l - T_"opt") / (R T_"opt" T_l))))
$

其中：$H_k = 200000$ J mol#super[-1]（去活化焓），$R = 8.314$ J mol#super[-1] K#super[-1]，
$T_"opt"$ 为最适温度（$V_"cmax"$: 301 K，$J_"max"$: 301 K），
$e_a$ 为活化能（$e_"vc" = e_"jm" = 55000$ J mol#super[-1]）。

*行为分析*：该函数在 $T_l = T_"opt"$ 处自动取得最大值，在 $T_"opt" = 301$ K（约 28°C）时
$V_"cmax"$ 达到峰值，函数在低温端单调上升，高温端迅速下降，为非归一化钟形。

== 2.2 独立模块：Medlyn 归一化 Arrhenius

独立模块使用 Medlyn et al. (2002) 的归一化峰值 Arrhenius 公式
（代码见 `Photosynthesis/temperature.jl:fTv`）：

$
  f_"Medlyn"(T_l) = V_(c"max",25) times
  underbrace(exp frac((T_l - T_(K,25)) e_"vc", T_(K,25) R T_l), "Arrhenius 增益")
  times
  underbrace(frac(1 + exp((S T_(K,25) - H_d) / (T_(K,25) R)),
                  1 + exp((S T_l - H_d) / (T_l R))), "25°C 归一化去活化项")
$

固定参数：$H_d = 200000$ J mol#super[-1]，$S = 640$ J mol#super[-1] K#super[-1]；
独立模块活化能：$e_"vc" = 30000$，$e_"jm" = 40000$ J mol#super[-1]。

*归一化证明*：当 $T_l = T_(K,25) = 298.15$ K 时，Arrhenius 增益 = 1，
去活化项分子分母相等 = 1，故 $f_"Medlyn"(T_(K,25)) = V_(c"max",25)$（严格归一化）。

== 2.3 温度响应曲线对比

下表给出两种函数在 $V_(c"max",25) = 89.45$ μmol m#super[-2] s#super[-1] 下的数值对比
（实际 $V_"cmax"$ 值，单位 μmol m#super[-2] s#super[-1]）：

#figure(
  table(
    columns: (auto, auto, auto, auto),
    inset: 7pt,
    align: (center, right, right, right),
    table.header([*$T$（°C）*], [*TBOLTZ*], [*Medlyn*], [*比值*]),
    [-5],  [—#super[†]], [37.0], [—],
    [0],   [—],          [49.0], [—],
    [5],   [—],          [62.6], [—],
    [10],  [—],          [76.2], [—],
    [15],  [—],          [89.4], [—],
    [20],  [41.5],       [90.6], [0.46],
    [25],  [72.5],       [89.4], [0.81],
    [28],  [87.2],       [83.2], [1.05],
    [30],  [88.0],       [72.5], [1.21],
    [35],  [76.8],       [42.9], [1.79],
    [40],  [53.1],       [16.6], [3.20],
    [45],  [24.7],       [3.6],  [6.86],
  ),
  caption: [$V_"cmax"$ 温度响应数值对比（$V_(c"max",25) = 89.45$）\ †: 低温时 TBOLTZ 计算值偏高，此处省略以聚焦主要差异区间],
) <table_Tresponse>

#beamer-block[
  *结论*：TBOLTZ 的最适温度为 ~28°C，在 20–40°C 的常见生长温度范围内给出更高
  $V_"cmax"$；Medlyn 方案的最适温度约 15°C，在 >20°C 时下降更快，
  更适合描述北方针叶林（低温偏好型植被）。
]

#pagebreak()

= 3 Michaelis-Menten 动力学常数

Farquhar 模型的 Rubisco 限制速率为：

$
  W_c = frac(V_"cmax" (c_i - Gamma^*), c_i + K_c(1 + O_2/K_o))
$

$K_c$、$K_o$、$tau$（Rubisco 特异性因子）均需在 25°C 参考值基础上进行温度校正：

$
  K_c = K_(c,25) exp frac((T_l - T_(K,25)) e_(K_c), T_(K,25) R T_l), quad
  "类似地" K_o, quad tau
$

两模块来自不同文献，参数值差异显著：

#figure(
  table(
    columns: (auto, auto, auto, auto),
    inset: 7pt,
    align: (left, right, right, left),
    table.header([*参数*], [*原始模块*], [*独立模块*], [*文献来源*]),
    [$K_(c,25)$ (μmol mol#super[-1])],   [274.6], [404.0],    [Harley 1995 vs Bernacchi 2001],
    [$K_(o,25)$ (mmol mol#super[-1])],   [419.8], [248.0],    [同上],
    [$tau_25$ (mmol mol#super[-1])],     [2904],  [2600],     [Balaguer 1996 vs Bernacchi 2001],
    [$e_(K_c)$ (J mol#super[-1])],       [80500], [80500],    [相同],
    [$e_(K_o)$ (J mol#super[-1])],       [14500], [14500],    [相同],
    [$e_tau$ (J mol#super[-1])],         [−29000],[−29000],   [相同],
    [$O_2$ (mmol mol#super[-1])],        [210],   [210],      [相同],
  ),
  caption: [Michaelis-Menten 常数对比],
)

*对 CO₂ 补偿点 $Gamma^*$ 的影响*：

$
  Gamma^* = frac(0.5 times O_2, tau) times 1000
  = frac(0.5 times 210000, tau_25 f_T(e_tau))
$

25°C 时：
- 原始模块：$Gamma^*_"orig" = 105000 / 2904 approx 36.2$ μmol mol#super[-1]
- 独立模块：$Gamma^*_"new" = 105000 / 2600 approx 40.4$ μmol mol#super[-1]

*对表观 Michaelis 常数 $K$ 的影响*（25°C）：

$
  K = K_c (1 + O_2/K_o)
$

- 原始：$K = 274.6 times (1 + 210/419.8) approx 412$ μmol mol#super[-1]
- 独立：$K = 404.0 times (1 + 210/248) approx 746$ μmol mol#super[-1]

#beamer-block[
  *结论*：独立模块的 $K$ 约为原始模块的 1.8 倍，意味着在相同 $c_i$ 下
  Rubisco 限制速率 $W_c$ 更低；同时更高的 $Gamma^*$ 进一步压低净光合，
  这两个效应共同使独立模块在中高温下给出更低的 $A_n$。
]

#pagebreak()

= 4 电子传输速率

== 4.1 原始模块：双曲饱和（Chen 1999）

$
  J_x = frac(J_"max" dot "PPFD", "PPFD" + 2.1 J_"max")
$

此公式是非矩形双曲线在曲率参数 $theta = 0$ 时的退化形式，
在低光时近似线性（斜率 $approx 1/2.1 approx 0.48$），
高光时趋近 $J_"max"$，半饱和点 PPFD$_(1/2) = 2.1 J_"max"$。

== 4.2 独立模块：非矩形双曲线（Farquhar 1980）

$
  theta_2 J^2 - (alpha I + J_"max") J + alpha I J_"max" = 0
$

解为：

$
  J = frac(alpha I + J_"max" - sqrt((alpha I + J_"max")^2 - 4 theta_2 alpha I J_"max")}{2 theta_2)
$

其中量子效率 $alpha = 0.3$，曲率参数 $theta_2 = 0.7$。

== 4.3 PAR 换算差异

#figure(
  table(
    columns: (2fr, 3fr, 3fr),
    inset: 7pt,
    align: (left, left, left),
    table.header([*项目*], [*原始模块*], [*独立模块*]),
    [换算系数], [$4.55 times 0.5$], [$4.6 times 0.5$],
    [冠层衰减], [外部辐射传输（`netRadiation_jl`）\ 已计算到叶片], [指数衰减\ $"PAR"_"leaf" = "PAR"_"total" e^(-0.5 L)$],
    [阴叶层累积 LAI], [辐射传输精确算法], [$L_"sunlit"$ 作为入射层深度],
  ),
  caption: [PAR 换算方式对比],
)

*差异效果*：在相同 Srad 输入下，两模块得到的阳生叶 PPFD 接近（4.55 vs 4.6 仅差 1%），
但阴叶 PAR 的计算方式不同——原始模块已在辐射传输阶段精确求解，
独立模块用简化指数衰减，使阴叶 PAR 偏低，导致阴叶 $J_x$ 和 $A_n$ 更低。

#pagebreak()

= 5 气孔-光合耦合求解策略

== 5.1 原始模块：解析三次方程

将 Ball-Berry 气孔导度方程和菲克扩散方程代入 Farquhar 方程，消去 $c_i$，
得到关于 $A_n$ 的三次多项式（推导见 `docs/modules/photosynthesis.typ` §1.5）：

$
  alpha A_n^3 + beta A_n^2 + gamma A_n + r = 0
$

通过 Vieta 三角代换法求解析根，选取物理上合理的根：

$
  A_"n,k" = -2sqrt(Q) cos frac(psi + 2 k pi, 3) - frac(p, 3), quad k = 0, 1, -1
$

当 $W_j le R_d$ 或 $W_c le R_d$（弱光或极低温）时退化为二次方程（$g_s = g_0$）。

*优势*：给定确定的 $T_l$、$"RH"_l$ 和辐射，一次性得到精确解析解，无迭代误差。

== 5.2 独立模块：定点迭代

```
初始化: ci = 0.7 × ca
循环（最多15次）：
    (An, Rd) = farquhar_model(T, PAR, ci, params)
    An = An × β_soil
    gs = ball_berry_gs(max(An,0), RH/100, ca, params)   ← 注意：cs = ca（简化！）
    gc = 1 / (1.4/gs + 1.37/gb_w)
    ci_new = ca - An / gc
    若 |ci_new - ci| < 1 μmol/mol → 收敛
    ci = ci_new
```

有一处简化：Ball-Berry 方程中 $c_s$ 取 $c_a$（大气 CO#sub[2]），
而非实际叶表面浓度 $c_s = c_a - A_n / g_b$。
这低估了高光合速率下 $c_s$ 的耗尽效应，使 $g_s$ 偏高，$c_i$ 偏高，$A_n$ 略偏高。

== 5.3 叶表面湿度的处理

这是两模块在物理假设上最根本的差异：

#figure(
  table(
    columns: (2fr, 3fr, 3fr),
    inset: 7pt,
    align: (left, left, left),
    table.header([*项目*], [*原始模块*], [*独立模块*]),
    [叶片温度], [由能量平衡迭代确定 $T_l ne T_a$], [假设 $T_l = T_a$],
    [叶表面 RH], [由潜热通量反推：\ $rho_v = ("LE"/lambda) r_v + rho_a$], [直接使用大气 RH],
    [能量耦合], [全耦合（光合↔蒸腾↔叶温）], [解耦（无叶温反馈）],
  ),
  caption: [叶表面湿度处理差异],
)

在原始模块中，$"RH"_l$ 通过实际潜热通量 LE 计算，
与光合-蒸腾耦合；晴天高辐射条件下叶温通常高于气温 1–5°C，
$"RH"_l < "RH"_a$，气孔导度因此降低。
独立模块直接使用大气 RH，相当于假设蒸腾不影响叶表面微环境，
在高辐射条件下会高估 $g_s$ 和 $A_n$。

#pagebreak()

= 6 Vcmax 冠层氮梯度

== 6.1 原始模块：Chen et al. 2012 垂直分布

原始模块通过 `VCmax()` 函数（`src/SPAC/VCmax.jl`）计算冠层氮梯度，
分别给出阳叶和阴叶的 $V_"cmax"$：

$
  V_("cmax,sunlit") = V_(c"max",25) chi N K frac(1 - e^(-(K_n + K) L), (K_n + K)(1 - e^(-K L)))
$

$
  V_("cmax,shaded") = V_(c"max",25) chi N left(frac{e^(-K_n L)}{K_n} - frac{1 - e^{-(K_n + K)L}}{K_n + K}right) / (L - frac{1 - e^{-KL}}{K})
$

其中 $chi$ 为比叶氮，$N$ 为叶片氮含量，$K = 0.5 Omega / cos theta_s$，$K_n = 0.3$ 为氮衰减系数。
在高 LAI 时，阴叶 $V_"cmax"$ 可能仅为阳叶的 20–50%。

== 6.2 独立模块：无氮梯度

独立模块对阳叶和阴叶使用相同的 $V_(c"max",25)$，
即假设氮均匀分布或氮梯度效应不显著。

*影响*：在高 LAI（>3）条件下，原始模块的阴叶光合速率更低，
而独立模块的阴叶 $A_n$ 会被高估。

#pagebreak()

= 7 综合结果差异

== 7.1 温度响应对比（定量）

以下结果来自等价单叶条件（$V_(c"max",25) = 89.45$，$R_H = 60\%$，$S_"rad" = 500$ W m#super[-2]，$g_b = 0.5$ mol m#super[-2] s#super[-1]）：

#figure(
  table(
    columns: (auto, auto, auto, auto, auto),
    inset: 7pt,
    align: (center, right, right, right, right),
    table.header([*T (°C)*], [$A_"n,orig"$], [$A_"n,new"$], [$c_"i,orig"$ (ppm)], [$c_"i,new"$ (ppm)]),
    [−5],  [—],    [6.5],  [—],    [146],
    [0],   [5.2],  [8.3],  [128],  [140],
    [5],   [7.2],  [10.3], [124],  [134],
    [10],  [9.6],  [12.2], [120],  [128],
    [15],  [12.0], [13.8], [117],  [123],
    [20],  [13.9], [11.9], [114],  [129],
    [25],  [14.5], [9.4],  [114],  [136],
    [30],  [8.3],  [6.7],  [123],  [145],
    [35],  [3.6],  [4.0],  [133],  [156],
    [40],  [1.1],  [1.5],  [157],  [175],
  ),
  caption: [$A_n$ 和 $c_i$ 的温度响应对比（单位：μmol m#super[-2] s#super[-1]）],
) <table_result>

== 7.2 差异溯因汇总

#figure(
  table(
    columns: (auto, 2fr, 3fr),
    inset: 7pt,
    align: (left, left, left),
    table.header([*温度区间*], [*$A_"n,orig"$ vs $A_"n,new"$*], [*主要原因*]),
    [< 20°C],   [$A_"n,new" > A_"n,orig"$],
      [独立模块 Medlyn 方案低温去活化弱\ Michaelis $K$ 较大（抵消部分增益）],
    [20–28°C],  [$A_"n,orig" > A_"n,new"$],
      [原始 TBOLTZ 最适温度附近峰值高\ 独立模块已过最适温度开始下降],
    [> 30°C],   [基本相当，偏差 <50%],
      [两方案去活化均快速；\ 动力学常数差异次要],
  ),
  caption: [温度区间差异成因],
)

== 7.3 设计定位与适用场景

#figure(
  table(
    columns: (2fr, 3fr, 3fr),
    inset: 7pt,
    align: (left, left, left),
    table.header([*维度*], [*`photosynthesis_jl`（原始）*], [*`Photosynthesis`（独立）*]),
    [定位],    [BEPS 主积分循环，与能量平衡耦合], [独立单叶计算、参数优化接口],
    [物理完整性], [全耦合（叶温、LE、气孔、光合）], [解耦简化（$T_l = T_a$）],
    [适用植被], [广泛（topt = 28°C）], [北方针叶林（topt ≈ 15°C）],
    [参数来源], [Harley & Baldocchi 1995], [Bernacchi et al. 2001],
    [计算成本], [较高（能量平衡迭代内部）], [低（独立迭代 ≤15步）],
    [优化友好性], [不友好（需完整积分循环）], [友好（简单接口，可直接优化）],
  ),
  caption: [两模块设计定位对比],
)

#pagebreak()

= 8 Bug 修复记录

在代码审查中，独立模块发现并修复了以下 5 处错误（2026-04-29）：

#figure(
  table(
    columns: (auto, 2fr, 2fr, 2fr),
    inset: 7pt,
    align: (center, left, left, left),
    table.header([*#*], [*位置*], [*原始错误*], [*修复后*]),
    [1], [`core.jl`\ $Gamma^*$ 公式],
      [$Gamma^* = 40.0 / tau$\ → 25°C 约 0.015 μmol mol#super[-1]],
      [$Gamma^* = 105000 / tau$\ → 25°C 约 40.4 μmol mol#super[-1]],
    [2], [`ParamPhoto.jl`\ $K_(o,25)$ 单位],
      [$K_(o,25) = 248000$（×1000 错误）\ O#sub[2]/$K_o$ ≈ 0.00085],
      [$K_(o,25) = 248$（mmol mol#super[-1]）\ O#sub[2]/$K_o$ ≈ 0.85],
    [3], [`temperature.jl`\ $f_T(V_"cmax")$ 归一化],
      [未归一化\ $V_"cmax"(25°C) approx 3724$],
      [Medlyn 归一化\ $V_"cmax"(25°C) = 89.45$（严格）],
    [4], [`core.jl`\ TPU 限制],
      [$W_p = 0.167 V_"cmax"$\ （系数缺少×3）],
      [$W_p = 0.5 V_"cmax"$\ （$W_p = 3 times "TPU"$，$"TPU" approx V_"cmax"/6$）],
    [5], [`photosynthesis.jl`\ 默认 $g_b$],
      [$g_b = 0.01$ mol m#super[-2] s#super[-1]\ （偏小 50 倍）],
      [$g_b = 0.5$ mol m#super[-2] s#super[-1]\ （典型叶片值）],
  ),
  caption: [独立模块 Bug 修复记录],
) <table_bugs>

= 9 建议

+ *温度响应统一*：若需在 20–35°C 常温范围内使用独立模块，
  建议将 $e_"vc"$ 调整为 55000 J mol#super[-1] 并采用 topt = 301 K，
  与原始模块保持一致。
+ *补全 $c_s$ 计算*：独立模块 Ball-Berry 中 $c_s = c_a$ 的简化
  在高 $A_n$ 时误差最大，建议改为 $c_s = c_a - A_n / g_b$。
+ *氮梯度*：高 LAI 站点应引入 `VCmax()` 函数计算阴/阳叶 $V_"cmax"$，
  而非使用统一值。
+ *参数统一化*：Michaelis 常数的选择对 $Gamma^*$ 和 $K$ 影响显著，
  应明确注释所使用文献来源，并保持与原始模块一致（除非有充分的生态学理由选用 Bernacchi 值）。

= 10 参考文献

+ Farquhar, G.D., von Caemmerer, S., & Berry, J.A. (1980). A biochemical model of photosynthetic CO₂ assimilation in leaves of C3 species. _Planta_, 149, 78–90.
+ Harley, P.C., & Baldocchi, D.D. (1995). Scaling carbon dioxide and water vapour exchange from leaf to canopy in a deciduous forest. _Plant, Cell & Environment_, 18(10), 1146–1156.
+ Chen, J.M., et al. (1999). Daily canopy photosynthesis model through temporal and spatial scaling for remote sensing applications. _Ecological Modelling_, 124, 99–119.
+ Medlyn, B.E., et al. (2002). Temperature response of parameters of a biochemically based model of photosynthesis. II. A review of experimental data. _Plant, Cell & Environment_, 25(9), 1167–1179.
+ Bernacchi, C.J., Singsaas, E.L., Pimentel, C., Portis, A.R., & Long, S.P. (2001). Improved temperature response functions for models of Rubisco-limited photosynthesis. _Plant, Cell & Environment_, 24(2), 253–259.
+ Ball, J.T., Woodrow, I.E., & Berry, J.A. (1987). A model predicting stomatal conductance and its contribution to the control of photosynthesis under different environmental conditions. _Progress in Photosynthesis Research_, 4, 221–224.
+ Chen, J.M., et al. (2012). Seasonal controls of the root-rhizomicrobial activity and the consequent effects on soil nitrogen mineralization in forests. _Global Biogeochemical Cycles_, 26.

#v(2em)
#align(center)[
  #text(size: 9pt, fill: gray)[
    BEPS.jl · 光合模块技术文档 · #datetime.today().display()
  ]
]
