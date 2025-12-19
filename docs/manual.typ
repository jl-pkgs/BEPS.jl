#import "@local/modern-cug-report:0.1.3": *
#show: doc => template(doc, footer: "BEPS.jl使用手册", header: "")

#import "@preview/codelst:2.0.2": sourcecode
// #set page(numbering: "1", number-align: center)
// #set text(font: "New Computer Modern", lang: "zh", size: 10.5pt)
// #set par(justify: true, first-line-indent: 2em)
// #show heading: set block(spacing: 1.5em)

#align(center)[
  #text(size: 24pt, weight: "bold")[BEPS.jl 模型文档]

  #v(1em)

  #text(size: 14pt)[Boreal Ecosystem Productivity Simulator]

  #v(0.5em)

  #text(size: 12pt)[Julia 实现版本]

  #v(2em)
]

#beamer-block[*目前这个手册是：Claude帮忙自动生成的，需谨慎使用。*]

#outline(indent: auto)

#v(1em)

#pagebreak()

= 1 模型概述


BEPS (Boreal Ecosystem Productivity Simulator) 是一个基于物理过程的陆地生态系统模型,用于模拟植被光合作用、蒸散发、土壤水热传输等关键生态水文过程。BEPS.jl 是该模型的 Julia 语言实现版本,具有高性能和良好的可扩展性。

== 1.1 主要特点

- *双层冠层结构*: 区分上层林冠 (overstorey) 和下层植被 (understorey),更精确地描述垂直结构对光、温、水的影响。
- *阴阳叶区分*: 考虑冠层内阳叶 (sunlit) 和阴叶 (shaded) 接受辐射的差异,分别计算光合作用和蒸腾。
- *完整的水热耦合*: 包含土壤水分运动(Richards方程)、热量传输(热传导方程)、冠层截留与蒸发等过程。
- *详细的辐射传输*: 区分直射和散射辐射,考虑冠层聚集指数对辐射传输的影响。
- *机理性光合模型*: 基于 Farquhar-von Caemmerer-Berry (FvCB) 模型,耦合气孔导度模型。
- *多层土壤模型*: 支持多层土壤水热模拟 (默认5层),能够模拟深层土壤过程。

== 1.2 模型结构

BEPS.jl 采用模块化设计,核心物理过程包括:
- *SPAC模块*: 模拟土壤-植物-大气连续体的能量和水分交换。
- *辐射模块*: 计算短波和长波辐射在冠层和地表的传输与分配。
- *光合模块*: 计算碳同化速率。
- *土壤模块*: 更新土壤温度、水分和积雪状态。

#pagebreak()

= 2 数据变量说明

模型主要处理以下几类关键数据:

== 2.1 气象驱动数据 (Meteorological Data)

模型运行需要逐小时或更高时间分辨率的气象驱动数据:

#table(
  columns: (1fr, 2fr, 1fr),
  inset: 8pt,
  align: horizon,
  table.header([*变量名*], [*描述*], [*单位*]),
  [temp], [空气温度], [℃],
  [rh], [相对湿度], [%],
  [rain], [降水量], [mm],
  [wind], [风速], [m/s],
  [Srad], [总太阳辐射], [W/$m^2$],
  [ca], [大气$"CO"_2$浓度], [ppm],
)

== 2.2 状态变量 (State Variables)

模型模拟过程中维护并更新的关键状态变量:

*植被状态*:
- *截留水/雪量*: 冠层叶片上截留的降雨和降雪量。
- *叶片温度*: 分层(上/下层)分叶(阴/阳叶)的叶片温度。

*土壤状态*:
- *地表温度*: 包括裸土温度、积雪表面温度等。
- *土壤温度*: 各层土壤的温度 $T_"soil"$.
- *土壤含水量*: 各层土壤的体积含水量 $theta$.
- *积雪状态*: 积雪深度 $z_"snow"$ 和积雪密度 $rho_"snow"$.
- *积水深度*: 地表积水深度 $z_"water"$.

== 2.3 冠层结构分层

模型将冠层垂直分为两层,水平分为两类叶片:

- *垂直分层*:
  - *Overstorey (上层)*: 乔木层或高大植被。
  - *Understorey (下层)*: 灌木、草本或苔藓层。
- *光照分叶*:
  - *Sunlit (阳叶)*: 接受直射和散射辐射。
  - *Shaded (阴叶)*: 仅接受散射辐射。

// #pagebreak()

= 3 模型参数

== 3.1 植被参数 (Vegetation Parameters)

不同植被类型(如针叶林、阔叶林、草地等)具有不同的生理生态参数:

#table(
  columns: (1.5fr, 3fr, 1fr),
  inset: 8pt,
  align: horizon,
  table.header([*参数符号*], [*描述*], [*典型值*]),
  [LAI_max], [最大叶面积指数], [2.0-6.0],
  [$alpha_"vis"$], [可见光反照率], [0.05-0.10],
  [$alpha_"nir"$], [近红外反照率], [0.20-0.40],
  [$Omega$], [聚集指数], [0.5-0.9],
  [$V_(c max 25)$], [25℃最大羧化速率], [30-100],
  [$g_1$], [气孔导度斜率参数], [5-12],
  [$z_"canopy"$], [冠层高度], [m],
)

== 3.2 土壤参数 (Soil Parameters)

土壤水力与热力性质取决于土壤质地(如砂土、壤土、粘土):

#table(
  columns: (1.5fr, 3fr, 1fr),
  inset: 8pt,
  align: horizon,
  table.header([*参数符号*], [*描述*], [*单位*]),
  [$theta_"sat"$], [饱和含水量], [$m^3/m^3$],
  [$theta_"vwp"$], [与凋萎点含水量], [$m^3/m^3$],
  [$K_"sat"$], [饱和导水率], [m/s],
  [$psi_"sat"$], [饱和水势], [m],
  [$b$], [Campbell参数], [-],
  [$kappa_"dry"$], [干土导热率], [W/m/K],
)

#pagebreak()

= 4 核心物理过程原理

== 4.1 太阳辐射传输算法

=== 4.1.1 太阳天顶角计算

太阳天顶角 ($theta_s$) 的余弦值决定了直射辐射的几何关系:

$ cos theta_s = cos delta cos phi cos H + sin delta sin phi $

其中:
- $delta$: 太阳赤纬角,通过儒略日计算
- $phi$: 纬度 (弧度)
- $H$: 太阳时角,与当地时间和经度相关

太阳赤纬角的计算采用傅里叶级数展开:

$
  delta = 0.006918 - 0.399912 cos((2 pi "DOY")/365) + 0.070257 sin((2 pi "DOY")/365) - \
  0.006758 cos((4 pi "DOY")/365) + 0.000907 sin((4 pi "DOY")/365)
$

=== 4.1.2 直射与散射辐射分离

总太阳辐射 ($R_s$) 首先被分解为直射 ($R_s^("dir")$) 和散射 ($R_s^("diff")$) 两个分量。模型采用经验公式基于云量比 ($r_c$) 来分离:

$ r_c = R_s / (1367 cos theta_s) $

$
  R_s^("diff") = cases(
    0.13 R_s & "if" r_c > 0.8,
    f(r_c) R_s & "otherwise"
  )
$

其中 $f(r_c)$ 是一个四次多项式拟合函数:

$ f(r_c) = 0.943 + 0.734 r_c - 4.9 r_c^2 + 1.796 r_c^3 + 2.058 r_c^4 $

直射辐射为残差: $R_s^("dir") = R_s - R_s^("diff")$

=== 4.1.3 冠层辐射传输

*直射辐射衰减*:

基于Beer-Lambert定律,考虑冠层聚集效应:

$ tau_("dir") = exp(-0.5 Omega "LAI" / cos theta_s) $

其中 $Omega$ 为冠层聚集指数 (0-1),反映叶片非随机分布。

*散射辐射衰减*:

散射辐射的传输使用等效天顶角:

$ cos Q = 0.537 + 0.025 "LAI" $

$ tau_("diff") = exp(-0.5 Omega "LAI" / cos Q) $

*阴阳叶LAI分配*:

阳叶LAI (接受直射光):
$ "LAI"_("sunlit") = (1 - exp(-0.5 Omega "LAI" / cos theta_s)) (2 cos theta_s) / Omega $

阴叶LAI (仅接受散射光):
$ "LAI"_("shaded") = "LAI" - "LAI"_("sunlit") $

=== 4.1.4 净辐射计算

*净短波辐射*计算考虑多次反射:

上层林冠:
$
  R_("ns")^o = R_s^("dir") [(1-alpha_o) - (1-alpha_u) tau_o^("dir")] + \
  R_s^("diff") [(1-alpha_o) - (1-alpha_u) tau_o^("diff")] + \
  0.21 Omega R_s^("dir") (1.1 - 0.1 "LAI"_o) exp(-cos theta_s)
$

*净长波辐射*:

基于Stefan-Boltzmann定律和多次反射:

$ R_("nl")^o = epsilon_o [R_l^("air") + R_l^u (1-tau_u^("diff")) + R_l^g tau_u^("diff")] - 2 R_l^o $

其中 $R_l = epsilon sigma T^4$ 为黑体辐射,总净辐射: $ R_n = R_("ns") + R_("nl") $

== 4.2 光合作用算法 (Farquhar模型)

=== 4.2.1 核心生化方程

光合速率受两个生化过程限制:

*1) Rubisco限制* ($W_c$,RuBP饱和):

$ W_c = (V_("cmax") (c_i - Gamma^*)) / (c_i + K) $

其中: $K = K_c (1 + O_2 / K_o)$ 为复合米氏常数。

*2) RuBP再生限制* ($W_j$,光限制):

$ W_j = (J (c_i - Gamma^*)) / (4 c_i + 8 Gamma^*) $

电子传递速率: $J = (J_("max") "PPFD") / ("PPFD" + 2.1 J_("max"))$

光合有效辐射: $"PPFD" = 4.55 times 0.5 times R_("ns")$ [$mu "mol" m^(-2) s^(-1)$]

*净光合速率*:

$ A_n = min(W_c, W_j) - R_d $

其中 $R_d = 0.015 V_("cmax")$ (光下呼吸减少40%)。

=== 4.2.2 参数的温度响应

关键参数随叶温变化:

*Arrhenius函数*:

$ f(T) = f_25 exp[(T - 25) E_a / (298.15 R T_K)] $

应用于: $K_c, K_o, tau, Gamma^*$

*峰值型函数*:

$
  f(T) = (f_25 (T_K / 298.15) exp[E_a (T_K - 298.15) / (298.15 R T_K)]) / \
  (1 + exp[(298.15 Delta S - H_d) / (298.15 R)])
$

应用于: $V_("cmax"), J_("max")$ (高温抑制)

=== 4.2.3 气孔导度模型

*Ball-Berry模型*:

$ g_s = g_0 + g_1 (A_n / c_s) "RH" beta_("soil") $

其中:
- $c_s = c_a - A_n / g_a$: 叶面$"CO"_2$浓度
- $beta_("soil")$: 土壤水分胁迫因子 (0-1),计算见4.4.4节
- $"RH"$: 叶面相对湿度

=== 4.2.4 光合-气孔耦合求解

这是一个非线性耦合系统。将Ball-Berry方程代入$"CO"_2$扩散方程:

$ c_i = c_a - A_n (1/g_a + 1/g_s) $

可得关于 $A_n$ 的*三次方程*:

$ alpha A_n^3 + beta A_n^2 + gamma A_n + theta = 0 $

*数值求解*采用Vieta替换法:

1. 令 $Q = (beta^2 - 3 gamma alpha) / (9 alpha^2)$
2. 计算 $U = (2 beta^3 - 9 alpha beta gamma + 27 alpha^2 theta) / (54 alpha^3)$
3. 三个实根: $A_("n,i") = -2sqrt(Q) cos((psi + 2pi i)/3) -beta/(3 alpha)$

其中 $psi = arccos(U / sqrt(Q^3))$,选择中间大小的正根。

== 4.3 能量平衡与蒸散发算法

=== 4.3.1 叶片能量平衡方程

每片叶片满足能量守恒:

$ R_n = H + "LE" $

*显热通量*:

$ H = rho_a c_p (T_("leaf") - T_a) g_h $

*潜热通量*:

$ "LE" = (rho_a c_p (e_s(T_("leaf")) - e_a)) / (gamma (1/g_s + 1/g_a)) $

其中 $gamma = c_p P / (0.622 lambda)$ 为干湿表常数。

=== 4.3.2 叶温迭代求解

由于 $e_s(T_("leaf"))$ 是温度的非线性函数(Clausius-Clapeyron方程),需要迭代:

1. 初始猜测: $T_("leaf")^((0)) = T_a - 0.5$℃
2. 计算饱和水汽压: $e_s(T_("leaf")^((k))) = 0.611 exp[(17.27 T_("leaf")) / (T_("leaf") + 237.3)]$ [kPa]
3. 计算能量通量: $H^((k))$, $"LE"^((k))$
4. 更新叶温: $T_("leaf")^((k+1)) = T_a + (R_n - "LE"^((k))) / (rho_a c_p g_h)$
5. 收敛判断: $|T_("leaf")^((k+1)) - T_("leaf")^((k))| < 0.02$℃

通常3-5次迭代即可收敛。

=== 4.3.3 蒸腾计算 (Penman-Monteith方程)

对于有气孔控制的蒸腾:

$
  "LE"_("trans") = (Delta R_n + rho_a c_p "VPD" g_a) / (Delta + gamma (1 + g_a / g_s))
$

其中:
- $Delta = (d e_s) / (d T) approx 4098 e_s / (T_a + 237.3)^2$: 斜率
- $"VPD" = e_s(T_a) - e_a$: 饱和水汽压差 [kPa]

蒸腾速率: $T = "LE"_("trans") / lambda$ [$"mm"/h$ 或 kg/$m^2$/s]

=== 4.3.4 截留蒸发与土壤蒸发

*冠层截留蒸发*:

湿叶表面蒸发(无气孔阻力):

$ E_("water") = (rho_a c_p (e_s(T_("leaf")) - e_a)) / (lambda gamma (1/g_a + 1/g_b + r_("min"))) $

其中 $r_("min") = 100$ s/m 为最小表面阻力。

*土壤蒸发*:

$ E_("soil") = (rho_a c_p (e_s(T_s) - e_a))/(lambda (r_a + r_s(theta))) $

土壤表面阻力: $r_s(theta) = r_("s,min") exp(beta (theta_("sat") - theta) / theta_("sat"))$

== 4.4 土壤水热传输算法

=== 4.4.1 Richards方程数值求解

*控制方程*:

$ (partial theta) / (partial t) = (partial) / (partial z) [K(theta) ((partial psi(theta)) / (partial z) + 1)] - S(z) $

*Campbell水力模型*:

$ psi(theta) = psi_("sat") (theta / theta_("sat"))^(-b) $

$ K(theta) = K_("sat") (theta / theta_("sat"))^(2b+3) $

*有限差分离散* (隐式格式):

$ (theta_i^(n+1) - theta_i^n) / (Delta t) = (F_(i-1/2)^(n+1) - F_(i+1/2)^(n+1)) / (Delta z_i) - S_i^(n+1) $

界面通量: $F_(i+1/2) = K_(i+1/2) ((psi_(i+1) - psi_i) / (Delta z_(i+1/2)) + 1)$

*界面导水率*采用加权几何平均:

$ K_(i+1/2) = (K_i psi_i + K_(i+1) psi_(i+1)) / (psi_i + psi_(i+1)) times (b_i + b_(i+1)) / (b_i + b_(i+1) + 6) $

*根系吸水分布*:

每层的吸水由经过水分胁迫修正的根系分布决定:

$ S_i = (T_o + T_u) d_i $

其中 $d_i$ 为第$i$层的有效根系吸水比例(见4.4.4节)。

=== 4.4.2 土壤热传导数值求解

*控制方程*:

$ C_s (partial T) / (partial t) = (partial) / (partial z) [kappa (partial T) / (partial z)] $

*体积热容*:

$ C_s = (1 - theta_("sat")) c_m rho_("soil") + theta c_w rho_w $

*土壤导热率* (Johansen模型):

干土: $kappa_("dry") = (0.135 rho_("soil") + 64.7) / (2700 - 0.947 rho_("soil"))$

饱和土: $kappa_("sat") = kappa_s^(1-theta_("sat")) kappa_w^(theta_("sat"))$

非饱和土: $kappa = (kappa_("sat") - kappa_("dry")) K_e + kappa_("dry")$

Kersten数: $K_e = log_10(S_r) + 1.0$ (当$S_r > 0.1$时)

其中 $S_r = (theta - theta_("res")) / (theta_("sat") - theta_("res"))$ 为相对饱和度。

*地表边界条件*:

$ T_g^(n+1) = (T_g^n C_s / Delta t + (R_n - "LE") + kappa_1 / Delta z_1 T_1^n) / (C_s / Delta t + kappa_1 / Delta z_1) $

=== 4.4.3 自适应时间步长

土壤水分计算采用自适应步长以平衡精度和效率:

$
  Delta t = cases(
    1.0 "s" & "if" max|F| > 10^(-5) "m/s",
    30.0 "s" & "if" max|F| > 10^(-6) "m/s",
    360.0 "s" & "otherwise"
  )
$

=== 4.4.4 土壤水分胁迫与根系分布

*根系垂直分布*:

BEPS采用指数衰减模型描述根系密度随深度的分布:

$ f_("root",i) = r^(100 z_(i-1)) - r^(100 z_i) $

其中:
- $z_i$ 为第$i$层土壤的累积深度 [m]
- $r$ 为根系衰减率参数 (典型值0.90-0.97)
- 因子100用于将深度从米转换为厘米

第一层: $f_("root",1) = 1 - r^(100 z_1)$

最后一层: $f_("root",N) = r^(100 z_(N-1))$ (吸收所有剩余根系)

*单层水分胁迫因子*:

每层土壤对根系吸水的限制通过水势计算:

$
  f_("psi",i) = cases(
    1 / (1 + ((psi_i - psi_("min")) / psi_("min"))^alpha) & "if" psi_i > psi_("min"),
    1.0 & "otherwise"
  )
$

其中:
- $psi_i = psi_("sat",i) (theta_i / theta_("sat",i))^(-b_i)$ 为该层水势
- $psi_("min")$ 为气孔关闭水势 (通常为10-33 m H₂O)
- $alpha$ 为水分胁迫指数 (0.4-1.5)

同时考虑温度对根系活性的影响:

$
  f_(T, i) = cases(
    1 - exp(-0.02 T_("soil", i)^2) & "if" T_("soil",i) > 0℃,
    0 & "otherwise"
  )
$

综合胁迫: $f_("stress",i) = f_("psi",i) times f_(T,i)$

*有效根系吸水分布*:

考虑土壤水分和温度胁迫后，根系吸水比例需重新归一化:

$ tilde(f)_i = f_("root",i) times f_("stress",i) $

$ d_i = tilde(f)_i / (sum_(j=1)^N tilde(f)_j) $

*整体水分胁迫因子*:

作用于气孔导度的整体土壤水分胁迫因子为:

$
  beta_("soil") = sum_(i=1)^N f_("stress",i) d_i = (sum_(i=1)^N f_("root",i) f_("stress",i)^2) / (sum_(j=1)^N f_("root",j) f_("stress",j))
$

该因子取值0.1-1.0,当 $beta_("soil") < 0.1$ 时强制设为0.1以避免气孔完全关闭。

*对光合和蒸腾的作用*:

- *光合作用*: $beta_("soil")$ 直接乘入Ball-Berry气孔导度方程,降低$g_s$,从而限制$c_i$和$A_n$
- *蒸腾速率*: 通过降低$g_s$,间接减小P-M方程中的$"LE"_("trans")$

== 4.5 积雪过程算法

=== 4.5.1 冠层截留雪

截获效率: $tau = 1 - exp(-"LAI" times Omega)$

截留雪量: $m_("snow")^(n+1) = m_("snow")^n + rho_("new") r_("snow") Delta t tau$

最大容量: $m_("max") = 0.1 "LAI"$ [kg/m²]

=== 4.5.2 雪密度演变

新雪密度: $rho_("new") = 67.9 + 51.3 exp(T_a / 2.6)$ [kg/m³]

压实(有降雪): $rho_("snow")^(n+1) = (rho_("snow")^n z^n + rho_("new") Delta z) / (z^n + Delta z)$

老化(无降雪): $rho_("snow")^(n+1) = (rho_("snow")^n - 250) exp(-0.001 Delta t / 3600) + 250$

=== 4.5.3 融雪能量平衡

融雪热量: $E_("melt") = m_("snow") c_("ice") T_("snow")$ [J/m²]

融雪量: $m_("melt") = min(E_("melt") / lambda_("fusion"), m_("snow"))$

其中 $lambda_("fusion") = 3.34 times 10^5$ J/kg

== 4.6 模型时间积分策略

BEPS采用*操作分裂法*:

1. *主时间步* (1小时): 更新气象驱动
2. *子时间步* (6分钟×10): 能量平衡迭代、蒸散发计算
3. *微观时间步* (自适应1-360秒): 土壤水分数值求解

这种多尺度设计在保证精度的同时最大化计算效率。

#pagebreak()

= 5 模型运行流程

BEPS模型按时间步长(通常为1小时)进行模拟,主要流程如下:

1. *初始化*: 读取配置参数、植被土壤参数,初始化各状态变量。
2. *时间循环*:
  - *读取气象数据*: 获取当前时刻的气温、降水、辐射等。
  - *降水分配*: 判断降雨/降雪,计算冠层截留量和穿透雨/雪量。
  - *能量平衡迭代*:
    这是一个核心的迭代过程,旨在求解冠层和地表的平衡温度。在迭代中同时计算:
    - 辐射分配
    - 光合作用速率
    - 气孔导度
    - 显热和潜热通量(蒸散发)
  - *状态更新*:
    - 计算实际蒸腾和土壤蒸发。
    - 更新冠层截留水量/雪量。
    - *土壤过程*: 求解热传导方程更新地温,求解Richards方程更新土壤水。
3. *输出结果*: 保存GPP、ET、土壤温湿度等结果。

#pagebreak()

= 6 使用说明

== 6.1 快速开始

BEPS.jl 提供了简洁的接口来运行模型。用户主要与 `besp_main` 函数交互。

```julia
using BEPS
using DataFrames

# 1. 准备气象驱动数据 (DataFrame格式)
# 必须包含: year, day, hour, tem, rh, rain, wind, Srad
drivers = DataFrame(
  year = fill(2020, 24),
  day  = fill(1, 24),
  hour = 1:24,
  tem  = 15.0 .+ randn(24),  # 气温
  rh   = 60.0 .+ rand(24),   # 湿度
  rain = zeros(24),          # 降水
  wind = fill(2.0, 24),      # 风速
  Srad = [0.0; range(0, 800, length=12); range(800, 0, length=11)] # 辐射
)

# 2. 准备LAI数据 (向量, 长度对应模拟天数)
lai = fill(3.0, 1) # 假设只有1天的LAI

# 3. 运行模型
# 可通过关键字参数指定经纬度、植被类型、土壤类型等
results, et, t_soil, theta_soil = besp_main(
  drivers, lai,
  lon = 120.0, lat = 30.0,
  VegType = 4,  # 植被类型代码 (如4=针叶林)
  SoilType = 8  # 土壤类型代码 (如8=壤土)
)
```

== 6.2 关键参数设置

在调用 `besp_main` 时,可以调整以下关键参数以适应不同站点的情况:

- `VegType`: 植被功能型 (Plant Functional Type) 代码。
  - 1: 常绿针叶林
  - 4: 落叶阔叶林
  - 5: 混交林
  - ... (详见Param模块)
- `SoilType`: 土壤质地代码。
  - 1-3: 砂土类
  - 4-6: 壤土类
  - 10-12: 粘土类
- `lat`, `lon`: 站点经纬度,用于计算太阳天顶角。
- `r_drainage`: 地表排水效率参数 (0-1),影响地表产流。
- `Tsoil0`, `θ0`: 土壤温度和湿度的初始条件。

== 6.3 结果变量说明

模型返回四个主要结果对象:

1. `df_out`: 包含主要的通量和状态结果。
  - `GPP`: 总初级生产力 ($"gC"/m^2$)
  - `Net_Rad`: 净辐射 ($W/m^2$)
  - `z_snow`: 积雪深度 (m)
2. `df_ET`: 详细的蒸散发分量。
  - `Trans_o/u`: 上/下层植被蒸腾 ($"mm"/h$)
  - `Eil_o/u`: 截留蒸发 ($"mm"/h$)
  - `Evap_soil`: 土壤蒸发 ($"mm"/h$)
3. `Tsoil`: 土壤温度矩阵 [时间 × 层数]。
4. `θ`: 土壤含水量矩阵 [时间 × 层数]。

#pagebreak()

= 7 参考文献

1. Chen, J.M., et al. (1999). Daily canopy photosynthesis model through temporal and spatial scaling for remote sensing applications. _Ecological Modelling_, 124(2-3), 99-119.
2. Liu, J., et al. (1997). A process-based boreal ecosystem productivity simulator using remote sensing inputs. _Remote Sensing of Environment_, 62(2), 158-175.
3. Farquhar, G.D., von Caemmerer, S., & Berry, J.A. (1980). A biochemical model of photosynthetic CO₂ assimilation in leaves of C3 species. _Planta_, 149(1), 78-90.
4. Ball, J.T., Woodrow, I.E., & Berry, J.A. (1987). A model predicting stomatal conductance and its contribution to the control of photosynthesis under different environmental conditions. _Progress in Photosynthesis Research_, 4, 221-224.
5. Campbell, G.S. (1974). A simple method for determining unsaturated conductivity from moisture retention data. _Soil Science_, 117(6), 311-314.

// #pagebreak()
#align(center)[
  // #text(size: 14pt, weight: "bold")[--- 文档结束 ---]
  #v(2em)

  #text(size: 10pt)[
    BEPS.jl 版本: v0.1.7 \
    文档生成日期: #datetime.today().display()
  ]
]
