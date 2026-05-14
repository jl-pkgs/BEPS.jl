#import "@preview/modern-cug-report:0.1.3": *
// #show: doc => template(doc, footer: "ModernBEPS.jl", header: "")

#show "{": "{"
#show "}": "}"


#pagebreak()
#let Delta(x) = $#sym.Delta #x$

// =
= 2 土壤温度

== 2.1 土壤温度模型结构

#h(2em)
ModernBEPS采用一维土壤热传导模型描述土壤温度变化。地表边界由净辐射、潜热通量、感热交换和积雪覆盖共同确定；土壤层内温度依据傅里叶热传导方程显式更新。模型默认划分为5层，层厚度为$Delta(z) = [0.05, 0.10, 0.20, 0.40, 1.25] "m"$。土壤温度模块的计算链条为：

$
  R_n, E, r_a, z_"snow" arrow.r "surface_temperature!"
  arrow.r G_1
  arrow.r "UpdateHeatFlux"
  arrow.r T_("soil",i), f_("ice",i)
$

#h(2em)
其中，`surface_temperature!`用于求解地表或雪表温度，并给出进入土壤系统的上边界热通量$G_1$；`UpdateHeatFlux`用于计算土壤层间热通量$G_i$，进而更新各层土壤温度。主要变量见表#[@table_module_ST]。

#figure(
  caption: [土壤温度模块主要变量。],
  table(
    columns: (2cm, 2.3cm, 3.2cm, 1fr, 2cm),
    align: (horizon, horizon, horizon, horizon, horizon),
    inset: 4pt,
    rows: 0.75cm,
    [*类别*], [*符号*], [*模型变量*], [*含义*], [*单位*],
    [State], [$T_("soil",i)^t$], [`Tsoil_p[i]`], [第 $i$ 层上一时间步土壤温度], [$degree "C"$],
    [State], [$T_("soil",i)^(t+1)$], [`Tsoil_c[i]`], [第 $i$ 层当前时间步土壤温度], [$degree "C"$],
    [State], [$f_("ice",i)$], [`ice_ratio[i]`], [第 $i$ 层土壤水冻结比例], [-],
    [State], [$G_i$], [`G[i]`], [第$i$个界面热通量，向下传输为正], [$"W m"^(-2)$],
    [State], [$C_(v,i)$], [`Cv[i]`], [第 $i$ 层体积热容], [$"J m"^(-3) "K"^(-1)$],
    [State], [$kappa_i$], [`κ[i]`], [第 $i$ 层土壤热导率], [$"W m"^(-1) "K"^(-1)$],
    // [Param], [$Delta(t)$], [`kstep`], [土壤温度更新步长，默认360 s], [$"s"$],
    [Param], [$rho_("soil",i)$], [`ρ_soil[i]`], [土壤容重], [$"kg m"^(-3)$],
    [Param], [$V_("SOM",i)$], [`V_SOM[i]`], [土壤有机质体积分数], [-],
    [Param], [$kappa_("dry",i)$], [`κ_dry[i]`], [干土热导率参数], [$"W m"^(-1) "K"^(-1)$],
    [Input], [$R_(n,g)$], [`radiation_g`], [地表净辐射], [$"W m"^(-2)$],
    [Input], [$g_(h,g)$], [`Gheat_g`], [地表空气动力学传热导度], [$"m s"^(-1)$],
    [Input], [$z_"snow"$], [`z_snow`], [积雪深度], [$"m"$],
  ),
) <table_module_ST>

#table-note[
  本章约定$G_i > 0$表示热量向下传输。`G[1]`为地表边界通量，`G[i]`（$i >= 2$）为第$i-1$层土壤底部界面通量。
]

=== 2.1.1 SnowLand 温度状态

#h(2em)
`SnowLand`用于保存地表与积雪相关温度状态。该结构并不只表示积雪温度，而是同时包含混合地表温度、雪面温度和雪下土壤表面温度。各变量含义见表#[@table_snowland]。

#figure(
  caption: [`SnowLand`温度状态变量。],
  table(
    columns: (2.7cm, 2.4cm, 1fr, 3.4cm),
    align: (horizon, horizon, horizon, horizon),
    rows: (0.9cm, 2.5cm, 2.5cm, 2.5cm, 1.3cm),
    fill: (x, y) => if calc.odd(y) { luma(97%) } else { none },
    inset: 4pt,
    [*模型变量*], [*符号*], [*含义*], [*情形*],

    [`T_snow0`],
    [$T_("snow",0)$],
    [*表层雪温度*。深雪情形下表示三层雪模型的最上层节点温度，并作为混合地表温度$T_("surf")$。],
    [薄雪、深雪],

    [`T_snow1`],
    [$T_("snow",1)$],
    [深雪三层模型中的*中层雪温度*。该层位于雪表层与雪底层之间，主要控制雪层内部热传导。],
    [深雪],

    [`T_snow2`],
    [$T_("snow",2)$],
    [深雪三层模型中的*底层雪温度*。该层连接雪层与雪下土壤表面，是计算雪-土界面热通量$G_("snow",2)$的上侧温度。],
    [深雪],

    [`T_mix0`], [$T_("mix",0)$], [*雪-土界面*加权温度。], [薄雪、深雪],

    [`T_surf`], [$T_("surf")$], [*地表与雪表*加权混合温度。], [所有情形],
  ),
) <table_snowland>

#table-note[
  `T_snow0`、`T_snow1`和`T_snow2`只在深雪条件下形成明确的三层雪温度廓线；在无雪或薄雪条件下，后两者主要用于保持状态变量完整，不代表独立雪层。
]

== 2.2 土壤热参数

=== 2.2.1 体积热容$C_(v,i)$

#h(2em)
每层土壤的体积热容由矿物颗粒、水/冰和有机质三部分组成。模型采用 Chen (2007) Eq. 18：

$
  C_(v,i)
  = 2.0 times 10^6 rho_("soil",i) / 2650
  + 10^6 theta_i [4.2(1 - f_("ice",i)) + 2.09 f_("ice",i)]
  + 2.5 times 10^6 V_("SOM",i)
$ <eq_Cv>

式中，第一项为矿物颗粒热容，第二项为土壤水和冰的热容，第三项为有机质热容。$theta_i$为第$i$层体积含水量，$f_("ice",i)$为冻结比例。

=== 2.2.2 土壤热导率$kappa_i$

#h(2em)
土壤热导率由干土参数、冰含量、水含量和相对饱和度共同决定：

$
  kappa_i =
  (kappa_("dry",i)^(1 - theta_("sat",i))
    k_i^(1.2 theta_i f_("ice",i))
    k_w^(theta_i (1 - f_("ice",i))) - 0.15)
  theta_i / theta_("sat",i) + 0.15
$ <eq_kappa_soil>

其中，$k_i = 2.1 "W m"^(-1) "K"^(-1)$为冰热导率，$k_w = 0.61 "W m"^(-1) "K"^(-1)$为水热导率。计算结果采用下限保护：

$ kappa_i arrow.r.long max(kappa_i, 0.15) $

该形式使热导率随含水量增加而增大；土壤冻结时，冰的较高热导率进一步增强传热。

=== 2.2.3 积雪热导率

#h(2em)
模型在求解积雪温度时用到积雪热导率$kappa_("dry_snow")$。
$kappa_("dry_snow")$由雪密度确定，对应代码变量`κ_dry_snow`。模型采用 Jordan (1991) 经验式：

$
  kappa_("dry_snow") = 0.021 + 4.2 times 10^(-4) rho_"snow"
  + 2.2 times 10^(-9) rho_"snow"^3
$ <eq_kappa_snow>

其中，$rho_"snow"$为积雪密度，单位为$"kg m"^(-3)$。

== 2.3 地表能量边界

=== 2.3.1 地表可用能量

#h(2em)
地表温度模块首先将净辐射和蒸发耗能合并为地表可用能量$G_g$：

$
  G_g = R_(n,g) - lambda_"snow" E_"snow" - lambda_w (E_"water" + E_"soil")
$ <eq_Gg>

其中，$E_"snow"$、$E_"water"$和$E_"soil"$分别为雪升华、地表水蒸发和土壤蒸发通量；$lambda_"snow" = 2.83 times 10^6 "J kg"^(-1)$，$lambda_w$为随气温变化的水汽化潜热。

#h(2em)
感热交换项可写为：

$ H_g = rho_a c_p (T_("surf") - T_a) / r_(a,g) $

其中，地表空气动力学阻力由传热导度给出$r_(a,g) = 1 / g_(h,g)$。

=== 2.3.2 表层控制体温度隐式求解

#h(2em)
令待求温度为$T^* = T^(t+1)$，上一时间步温度为$T^t$。这里的$T^*$不是固定意义上的“两层交界处温度”，而是`solve_imp`当前要求解的表层控制体温度。根据调用位置不同，它可以表示裸土地表温度、雪表层温度，或与雪土界面相关的表层温度。该温度变化由三类通量共同决定：一是地表可用能量输入，二是与上边界空气或冠层之间的感热交换，三是向下方介质的热传导。通过求解以下形式的能量平衡方程，可得到其解：

$ C_v Delta(z) pdv(T, t) = G_"net" - H - G_"down" $
其中，$G_"net" = R_n - "LE"$ 为净辐射及其他源项，$G_"down"$ 为向下传导热通量。

// 控制体的储热项为$C_v (Delta z) / Delta(t) (T^* - T^t)$。其中，$C_v Delta z$为单位面积热容。

#h(2em)
等号右侧的三个项分别为：地表可用能量输入、上边界感热损失和向下传导损失。为避免混淆，下文用$Delta z_M$表示$T^*$所代表控制体的储热厚度，用$z_c$表示`solve_imp`中传入的交换距离尺度。二者有时相等，但在当前代码中并不总是同一个量。

三个通量项可写作：
$
  "地表可用能量:" G = G_g dot z_"rad" / z_c \
  "上边界感热损失为:" H = rho_a c_p (T^* - T_"up") / r_a \
  "向下传导损失为：" G_"down" = kappa_"bot" (T^* - T_"bot") / d_c = eta_c kappa_"bot" (T^* - T_"bot") / z_c, (eta_c = z_c / d_c)
$

其中，$d_c$为$T^*$到下边界温度$T_"bot"$之间的实际传导距离，$eta_c$为等效传导距离修正系数，对应模型代码中的`η_c`参数。薄雪情形的裸土地表温度求解使用`η_c=2.0`，表示地表到第一层土壤中心的传导距离取半层厚度。地表可用能量$G_g$不一定作用于整个交换距离$z_c$，因此引入作用深度$z_"rad"/z_c$修正因子。

#h(2em)
将上述各项带入，可得：

$
  Delta(M) (T^* - T^t)
  = G_g z_"rad" / z_c
  - rho_a c_p (T^* - T_"up") / r_a
  - eta_c kappa_"bot" (T^* - T_"bot") / z_c, quad Delta(M) = C_v Delta(z_M) / Delta(t)
$ <eq_imp_balance>

#h(2em)
$Delta(z_M)$不是“交界处的厚度”，而是$T^*$代表的控制体参与储热的厚度；$z_c$对应`solve_imp`参数`z`，表示该控制体与下方温度$T_"bot"$进行热交换时采用的距离尺度。因此，式#[@eq_imp_balance]中的$z_c$不应统一改写为$Delta(z_M)$；只有当储热厚度与交换距离相同时，才有$z_c = Delta(z_M)$。

#let dz = "dz"
// #beamer-block[
//   *示例：深雪雪表层温度*。在深雪情形中，雪表层温度`T_snow0`由以下代码求解：

//   ```julia
//   ΔM_snow = cp_ice * ρ_snow * dz_snow_s1 / Δt
//   T_snow0 = solve_imp(..., ra_g, dz_snow_s12, Gg, ρCp, κ_dry_snow;
//                        z_rad=dz_snow_s1)
//   ```

//   其中，$dz_"snow,s1" = 0.02 "m"$，$dz_"snow,s12" = dz_"snow,s1" + dz_"snow,s2" = 0.04 "m"$。因此：

//   $
//     Delta(z_M) = dz_"snow,s1" = 0.02 "m", quad
//     z_c = dz_"snow,s12" = 0.04 "m", quad
//     z_"rad" = dz_"snow,s1" = 0.02 "m"
//   $

//   这表示：雪表层只有$0.02 "m"$厚的雪参与储热，但其与下方雪层交换热量时采用$0.04 "m"$的距离尺度；同时，地表可用能量只作用在雪表层$0.02 "m"$厚度上。
// ]

#h(2em)
为得到$T^*$的显式表达式，将式#[@eq_imp_balance]两边同乘$r_a z_c$，消去分母：

$
  Delta(M) r_a z_c (T^* - T^t)
  = G_g r_a z_"rad"
  - rho_a c_p z_c (T^* - T_"up")
  - eta_c r_a kappa_"bot" (T^* - T_"bot")
$ <eq_imp_balance_scaled>

令$I = Delta(M) r_a z_c$，并展开各项：

$
  I T^* - I T^t =
  G_g r_a z_"rad" - rho_a c_p z_c T^* + rho_a c_p z_c T_"up"
  - eta_c r_a kappa_"bot" T^* + eta_c r_a kappa_"bot" T_"bot"
$

将所有含$T^*$的项移到左侧，其余已知项移到右侧：

$
  T^* (I + rho_a c_p z_c + eta_c r_a kappa_"bot")
  = I T^t
  + G_g r_a z_"rad"
  + rho_a c_p z_c T_"up"
  + eta_c r_a kappa_"bot" T_"bot"
$ <eq_imp_balance_collect>

因此，隐式温度解为：

$
  T^* =
  (T^t I + G_g r_a z_"rad" + rho_a c_p T_"up" z_c + eta_c r_a kappa_"bot" T_"bot")
  /
  (rho_a c_p z_c + eta_c r_a kappa_"bot" + I)
$ <eq_solve_imp>

其中，$T_"up"$为上边界空气或下层冠层温度，$T_"bot"$为下方介质温度，$kappa_"bot"$为下方介质热导率。式#[@eq_solve_imp]对应`solve_imp`函数。若给定参考温度$mu$，求解结果将截断至$[mu - 25, mu + 25]$，以抑制地表温度的非物理跳变。

== 2.4 积雪控制下的地表温度

#h(2em)
模型按积雪深度$z_"snow"$将地表边界划分为三种情形，并据此确定混合地表温度$T_("surf")$、雪-土界面温度$T_("mix",0)$和进入土壤的热通量$G_1$。

=== 2.4.1 无雪或极薄雪

#h(2em)
当$z_"snow" <= 0.02 "m"$时，地表视为无雪或极薄雪覆盖。此时地表温度、土壤表面温度和积雪相关温度取同一值：

$
  T_("surf") = T_("mix",0) = T_("soil",0) = T_("snow",0)
$

地表温度由式#[@eq_solve_imp]求解，下边界取第一层土壤温度$T_("soil",1)^t$。进入土壤的热通量为：

$
  G_1 = kappa_1 (T_("surf") - T_("soil",1)^t) / (Delta(z_1) \/ 2 )
$ <eq_G_surface_nosnow>

并截断至$[-100, 100] " W m"^(-2)$。

=== 2.4.2 薄雪或斑块雪

#h(2em)
当$0.02 "m" < z_"snow" <= 0.05 "m"$时，模型采用雪-裸土双板块方法：裸土部分求解裸土地表温度$T_("soil",0)$，积雪覆盖部分求解雪表温度$T_("snow",0)$。

#h(2em)
雪覆盖区域的雪-土界面温度$T_"int"$由土壤传导、积雪传导和界面热惯性共同决定：

#h(2em)
为说明式#[@eq_Tint]的来源，取垂直坐标$z$向下为正，并将雪-土界面视为带有等效储热项的薄层控制体。设界面温度为$T_"int"$，上一时间步界面温度为$T_("mix",0)^t$。雪表层温度为$T_("snow",0)$，位于界面上方；第一层土壤温度为$T_("soil",1)^t$，位于界面下方。

#h(2em)
傅里叶定律在向下为正的坐标下写为：

$
  G = - kappa (dif T) / (dif z)
$ <eq_fourier_down>

界面上侧雪层通量取“流向界面”为正。雪表层到界面的温度梯度近似为$(T_"int" - T_("snow",0)) / z_"snow"$，因此从雪层流入界面的通量为：

#let snow2int = $"snow" -> "int"$
#let int2soil = $"int" -> "soil"$

$
  G_snow2int
  = - kappa_("dry_snow") (T_"int" - T_("snow",0)) / z_"snow"
  = kappa_("dry_snow") (T_("snow",0) - T_"int") / z_"snow"
$ <eq_G_snow_int>

热通量向下为正。界面到第一层土壤节点的传导距离为$Delta(z_("soil"))$，因此从界面流出的通量为：

$
  G_int2soil
  = - kappa_1 (T_("soil",1)^t - T_"int") / Delta(z_("soil"))
  = kappa_1 (T_"int" - T_("soil",1)^t) / Delta(z_("soil"))
$ <eq_G_int_soil>

界面储热等于流入通量减去流出通量：

$
  Delta(M_1) (T_"int" - T_("mix",0)^t) =
  G_snow2int - G_int2soil
$

代入式#[@eq_G_snow_int]和式#[@eq_G_int_soil]：

$
  Delta(M_1) (T_"int" - T_("mix",0)^t) =
  kappa_("dry_snow") (T_("snow",0) - T_"int") / z_"snow"
  + kappa_1 (T_("soil",1)^t - T_"int") / Delta(z_("soil"))
$ <eq_Tint_balance>

其中，$Delta(M_1) = C_(v,1) dot 0.02 / Delta(t)$。$0.02 "m"$为人为引入的等效界面储热层；它使界面温度具有时间记忆，但并非真实几何界面厚度。

#h(2em)
将式#[@eq_Tint_balance]中所有含$T_"int"$的项移至左侧，可得：

$
  T_"int" =
  (kappa_1 T_("soil",1)^t / Delta(z_("soil"))
  + kappa_("dry_snow") T_("snow",0) / z_"snow"
  + Delta(M_1) T_("mix",0)^t)
  /
  (kappa_1 / Delta(z_("soil"))
  + kappa_("dry_snow") / z_"snow"
  + Delta(M_1))
$ <eq_Tint>

其中，$Delta(z_("soil")) = Delta(z_1) / 2$，$Delta(M_1) = C_(v,1) 0.02 / Delta(t)$。面积平均土壤表面温度为：

$
  T_("mix",0) = (1 - f_"snow") T_("soil",0) + f_"snow" T_"int"
$ <eq_Tmix_patchy>

上边界混合辐射温度为：

$
  T_("surf") = (1 - f_"snow") T_("soil",0) + f_"snow" T_("snow",0)
$ <eq_Tsurf_patchy>

进入土壤的热通量按积雪覆盖区和裸土区面积加权：

$
  G_1 =
  f_"snow" kappa_("dry_snow") (T_("snow",0) - T_("soil",1)^t) / (z_"snow" + Delta(z_("soil")))
  + (1 - f_"snow") kappa_1 (T_("soil",0) - T_("soil",1)^t) / Delta(z_("soil"))
$ <eq_G_surface_patchy>

其中，$f_"snow"$为地表积雪覆盖比例。

=== 2.4.3 深雪

#h(2em)
当$z_"snow" > 0.05 "m"$时，积雪主导地表边界。模型采用三层积雪温度节点（表#[@table_deepsnow_layers]）。

#[
  #show table.cell: it => {
    set par(leading: 0.8em, spacing: 0.5em)
    set text(size: 11.5pt)
    set math.text(size: 11pt)
    it
  }

  #figure(
    caption: [
      三层积雪相关变量。
    ],
    table(
      columns: (2.1cm, 2.4cm, 3.4cm, 3.2cm, 1fr),
      align: (horizon, horizon, horizon, horizon, horizon),
      rows: (0.95cm, 1.3cm),
      fill: (x, y) => if calc.odd(y) { luma(97%) } else { none },
      // inset: 4pt,
      [*层位*], [*温度变量*], [*厚度*], [*储热项*], [*主要热通量*],
      [雪表层],
      [$T_("snow",0)$],
      [$d_("snow0") = 0.02 "m"$],
      [$c_("ice") rho_"snow" d_("snow0")$],
      [受$G_g$、感热交换影响；$G_"snow,0"$表层向中层传导],

      [雪中层],
      [$T_("snow",1)$],
      [$d_("snow1") = 0.02 "m"$],
      [$c_("ice") rho_"snow" d_("snow1")$],
      [流入$G_"snow,0"$，流出$G_"snow,1"$],

      [雪底层],
      [$T_("snow",2)$],
      [$d_("snow2") = z_"snow" - 0.04 "m"$],
      [$c_("ice") rho_"snow" d_("snow2")$],
      [流入$G_"snow,1"$，流出$G_"snow,2"$],

      [雪下土壤表面],
      [$T_("mix",0)$],
      [$d_("soil,0") = 0.02 "m"$],
      [$C_(v,1) d_("soil,0")$],
      [流入$G_"snow,2"$，流出$G_1^t$],
    ),
  ) <table_deepsnow_layers>
]

雪表温度仍由式#[@eq_solve_imp]隐式求解；中层、底层和雪下土壤表面温度采用显式欧拉格式更新：

$
  T^(t+1) = T^t + (F_"in" - F_"out") Delta(t) / C
$ <eq_step_exp>

该情形下主要热通量为：

$
  G_"snow,0" = kappa_("dry_snow") (T_("snow",0) - T_("snow",1)^t) / 0.04
$

$
  G_"snow,1" = kappa_("dry_snow") (T_("snow",1)^t - T_("snow",2)^t) / (z_"snow" - 0.02)
$

$
  G_"snow,2" =
  (T_("snow",2)^t - T_("mix",0)^t)
  /
  [0.5 (z_"snow" - 0.04) / kappa_("dry_snow") + 0.02 / kappa_1]
$ <eq_G_snow2>

雪下土壤表面温度更新为：

$
  T_("mix",0)^(t+1)
  = T_("mix",0)^t
  + (G_"snow,2" - G_1^t) Delta(t) / (C_(v,1) dot 0.02)
$ <eq_Tmix_deepsnow>

其中，$G_1^t$为上一时间步进入第一层土壤的热通量。深雪条件下，混合地表温度取雪表温度：

$ T_("surf") = T_("snow",0) $

// #beamer-block[
//   *深雪节点距离说明*：本节中的雪层传导距离与 `surface_temperature.jl` 当前实现保持一致。`docs/ch02_SoilTemperature/TS_Case03.md` 给出了基于节点中心网格的距离审查推导，可作为后续修订深雪传导距离的依据。
// ]

== 2.5 土壤层内热传导

=== 2.5.1 层间热通量

#h(2em)
确定地表边界热通量$G_1$后，土壤内部热通量按傅里叶定律计算。第$i$层与第$i+1$层之间的界面通量为：

$
  G_(i+1) =
  2 (T_("soil",i)^t - T_("soil",i+1)^t)
  /
  (Delta(z_i) / kappa_i + Delta(z_(i+1)) / kappa_(i+1)),
  quad i = 1, dots, N - 1
$ <eq_G_internal>

该式等价于两层半层热阻串联后的界面传导。若上层温度高于下层，则$G_(i+1) > 0$，表示热量向下传输。

=== 2.5.2 底边界通量

#h(2em)
最底层连接深层恒温边界，深层温度以年平均气温$T_("air,ann")$近似。底边界通量为：

$
  G_(N+1) =
  kappa_N (T_("soil",N)^t - T_("air,ann"))
  /
  (D_f + Delta(z_N) / 2)
$ <eq_G_bottom>

其中，$D_f$对应`DEPTH_F`。土壤内部界面通量计算后截断至$[-200, 200] "W m"^(-2)$。

=== 2.5.3 温度更新

#h(2em)
对第$i$层土壤，热量守恒方程为：

$
  C_(v,i) Delta(z_i) (T_("soil",i)^(t+1) - T_("soil",i)^t) / Delta(t) = G_i - G_(i+1) + S_i^H
$ <eq_T_balance>

当前实现取内部热源$S_i^H = 0$，因此温度显式更新为：

$
  T_("soil",i)^(t+1) = T_("soil",i)^t
  + (G_i - G_(i+1)) Delta(t) / (C_(v,i) Delta(z_i))
$ <eq_T_update>

更新后，土壤温度截断至：

$ T_("soil",i)^(t+1) arrow.r.long "clamp"(T_("soil",i)^(t+1), -50, 50) $

若提供观测土壤温度并设置`fix_Tsoil=true`，模型将跳过式#[@eq_T_update]的温度更新，但仍更新冻融比例。

== 2.6 冻融状态更新

#h(2em)
土壤冻融状态以冻结比例$f_("ice",i)$表示。若第$i$层土壤温度由非负跨越至负值，且该层尚未完全冻结，则低于$0 degree "C"$的冷却能量用于冻结液态水：

$
  Q_f = (0 - T_("soil",i)^(t+1)) C_(v,i) Delta(z_i)
$

$
  f_("ice",i) arrow.r.long
  min(1, f_("ice",i) + Q_f / (L_f rho_w theta_i Delta(z_i)))
$ <eq_freeze>

随后令$T_("soil",i)^(t+1) = 0 degree "C"$。

#h(2em)
若温度由非正跨越至正值，且该层存在冰，则升温能量用于融冰：

$
  Q_m = (T_("soil",i)^(t+1) - 0) C_(v,i) Delta(z_i)
$

$
  f_("ice",i) arrow.r.long
  max(0, f_("ice",i) - Q_m / (L_f rho_w theta_i Delta(z_i)))
$ <eq_melt>

同样令$T_("soil",i)^(t+1) = 0 degree "C"$。模型中$L_f = 3.34 times 10^5 "J kg"^(-1)$；$rho_w = 1000 "kg m"^(-3)$在冻结比例换算中以代码形式体现。

#h(2em)
考虑水分更新引起的含水量变化，冻结比例进一步按上一时间步与当前含水量缩放：

$
  f_("ice",i) arrow.r.long
  min(1, f_("ice",i) theta_i^t / theta_i^(t+1))
$ <eq_ice_adjust>

该处理使冻结水量与土壤水储量变化保持一致，并将冻结比例限制在$[0, 1]$范围内。

== 2.7 计算流程小结

#h(2em)
每个次小时步长内，土壤温度模块按以下顺序计算：

$
  "UpdateThermal"_kappa
  arrow.r "UpdateThermal"_C_v
  arrow.r "surface_temperature!" \
  arrow.r G_1
  arrow.r "UpdateHeatFlux"
  arrow.r "Update_ice_ratio"
$

#h(2em)
模型首先根据当前含水量和冻结比例更新热导率与热容；随后由地表能量平衡和积雪状态计算上边界热通量；再通过层间热传导更新土壤温度；最后依据温度是否跨越$0 degree "C"$更新冻结比例。土壤温度还影响土壤水分模块中的冻结因子和根系吸水温度胁迫，由此形成热-水耦合。

// function surface_temperature_jl!(
// Rn_g::FT, T_air::FT, Tc_u::FT, RH::FT,
// z_snow::FT, z_water::FT, ρ_snow::FT, perc_snow_g::FT,
// z_soil1::FT, κ_soil1::FT, Cv_soil1::FT, Cv_soil0::FT, Gheat_g::FT,
// E_soil::FT, E_water_g::FT, E_snow_g::FT,
// G_soil1::FT, T_soil1_last::FT, T_soil0_last::FT,
// last::SnowLand{FT}, current::SnowLand{FT};
// kstep=360.0
// )

#pagebreak()
