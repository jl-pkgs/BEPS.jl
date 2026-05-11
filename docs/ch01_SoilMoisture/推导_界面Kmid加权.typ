#import "@preview/modern-cug-report:0.1.3": *
#show: doc => template(doc, footer: "CUG水文气象学2025", header: "")

= Campbell 型界面导水率：异质层界面的 $b$ 如何取值

本文只讨论物理和数值离散本身，不以 BEPS 既有方案为准。结论先写在前面：

- 对同一种土壤，Campbell 积分会给出 $b / (b + 3)$。
- 对相邻两层不同土壤，$b_i$ 和 $b_(i+1)$ 分别属于各自土层的保水曲线；界面处没有一个天然的、静态的 $b_"eff"$。
- 若做有限体积离散，界面两侧半层应分别使用自己的 $b_i$ 和 $b_(i+1)$，再通过通量连续或阻力串联得到界面通量。
- 因而 $b_"eff" = (b_i + b_(i+1)) / 2$ 只是经验平均，不是由 Campbell 积分推出的结果。

== 1. 层中心量与界面量

第 $i$ 层中心的 Campbell 关系为：

$
  psi_i = psi_("sat",i)
  (theta_i / theta_("sat",i))^(-b_i)
$

$
  K_i =
  K_("sat",i)
  (theta_i / theta_("sat",i))^(2 b_i + 3)
$

这里 $b_i$ 是第 $i$ 层材料的孔径分布参数。它描述的是这一层内部 $theta$、$psi$、$K$ 之间的函数关系，不是两个层中心之间可以直接插值的状态变量。

界面导水率 $K_(i+1/2)$ 是用于层间通量的数值等效量：

$
  q_(i+1/2) = K_(i+1/2) (
    frac(psi_(i+1) - psi_i, Delta z_(i+1/2)) + 1
  )
$ <eq_q>

这里的括号内有两项：

$
  q_(i+1/2)
  =
  underbrace(
    K_(i+1/2) frac(psi_(i+1) - psi_i, Delta z_(i+1/2)),
    "毛管势梯度项"
  )
  +
  underbrace(K_(i+1/2), "重力项")
$

因此，$+1$ 不是被吃掉了，而是对应重力项 $q_g = K$。下面的 $integral K dif psi$ 只能推出毛管项所需的界面导水率；它不能直接推出重力项。先只看毛管势梯度项：

$ q_psi = K(psi) frac(dif psi, dif z) $

在第 $i$ 层中心到第 $i+1$ 层中心之间，若用一个常数界面导水率表示同一段通量，则有限差分形式为：

$
  q_psi
  =
  K_(i+1/2)
  frac(psi_(i+1) - psi_i, Delta z_(i+1/2))
$

两边乘以 $Delta z_(i+1/2)$，得到：

$
  q_psi Delta z_(i+1/2)
  =
  K_(i+1/2) (psi_(i+1) - psi_i)
$

另一方面，由连续形式有：

$
  q_psi dif z = K(psi) dif psi
$

沿 $psi_i$ 到 $psi_(i+1)$ 积分：

$
  q_psi Delta z_(i+1/2)
  =
  integral_(psi_i)^(psi_(i+1)) K(psi) dif psi
$

把这两个表达式相等，便得到界面导水率的理论定义：

$
  K_(i+1/2)
  =
  frac(
    integral_(psi_i)^(psi_(i+1)) K(psi) dif psi,
    psi_(i+1) - psi_i
  )
$ <eq_k_interface>

这个定义保证毛管项满足：

$
  K_(i+1/2) (psi_(i+1) - psi_i)
  =
  integral_(psi_i)^(psi_(i+1)) K(psi) dif psi
$

若把重力项也写得更严格，应允许毛管项和重力项使用不同的界面平均：

$
  q_(i+1/2)
  =
  K_(psi,i+1/2)
  frac(psi_(i+1) - psi_i, Delta z_(i+1/2))
  +
  K_(g,i+1/2)
$

其中：

$
  K_(psi,i+1/2)
  =
  frac(
    integral_(psi_i)^(psi_(i+1)) K(psi) dif psi,
    psi_(i+1) - psi_i
  )
$

而 $K_(g,i+1/2)$ 应来自沿空间路径的平均：

$
  K_(g,i+1/2)
  =
  frac(1, Delta z_(i+1/2))
  integral_(z_i)^(z_(i+1)) K(z) dif z
$

把同一个 $K_(i+1/2)$ 同时放进毛管项和重力项，是一种简化离散；不是 $integral K dif psi$ 本身要求的结果。

其中 $K_i$ 是层中心导水率，$K_(i+1/2)$ 是界面等效导水率，二者不能混淆。

== 2. 同质层内的 Campbell 积分

先看同一种土壤，即 $b$、$psi_("sat")$、$K_("sat")$ 在这一段内不变。令：

$
  x = theta / theta_("sat")
$

Campbell 关系为：

$
  psi = psi_("sat") x^(-b),
  quad
  K = K_("sat") x^(2b + 3)
$

对 $psi$ 求微分：

$
  dif psi = - b psi_("sat") x^(-b - 1) dif x
$

于是：

$
  K dif psi
  =
  K_("sat") x^(2b + 3)
  (- b psi_("sat") x^(-b - 1)) dif x
  =
  - b K_("sat") psi_("sat") x^(b + 2) dif x
$

积分得到：

$
  integral K dif psi = - frac(b, b + 3) K_("sat") psi_("sat") x^(b + 3)
$

又因为：

$
  K psi = K_("sat") x^(2b + 3) dot psi_("sat") x^(-b)
  = K_("sat") psi_("sat") x^(b + 3)
$

所以：
$ integral K dif psi = - frac(b, b + 3) K psi $

由此可以定义同质层内的理论界面导水率。设第 $i$ 层中心到第 $i+1$ 层中心之间属于同一种土壤，且二者都落在同一条 Campbell 曲线上。令：

$ F(psi) = integral K(psi) dif psi = - frac(b, b + 3) K(psi) psi $

则界面导水率应取沿 $psi$ 路径的积分平均：

$
  K_(i+1/2)^"int" = frac(
    integral_(psi_i)^(psi_(i+1)) K(psi) dif psi,
    psi_(i+1) - psi_i
  )
  = frac(F(psi_(i+1)) - F(psi_i), psi_(i+1) - psi_i)
$

代入 $F(psi)$ 得：

$
  K_(i+1/2)^"int"
  =
  - frac(b, b + 3)
  frac(
    K_(i+1) psi_(i+1) - K_i psi_i,
    psi_(i+1) - psi_i
  )
$

当 $psi_(i+1) arrow.r psi_i$ 时，上式按极限退化为 $K_i$。因此，$b / (b + 3)$ 只在“沿同一条 Campbell 曲线积分”时成立；这个 $K_(i+1/2)^"int"$ 也只是在同质段内成立的理论界面导水率。

== 3. 异质界面不能直接取一个 $b$

若第 $i$ 层和第 $i+1$ 层土壤不同，则两侧有两条不同的水力曲线：

$
  K_i(psi) =
  K_("sat",i)
  (psi / psi_("sat",i))^(-(2 b_i + 3) / b_i)
$

$
  K_(i+1)(psi) =
  K_("sat",i+1)
  (psi / psi_("sat",i+1))^(-(2 b_(i+1) + 3) / b_(i+1))
$

界面处应连续的是压力水头和通量，而不是 $theta$，也不是 $b$。因此在材料界面上：

- $b_i$ 控制第 $i$ 层中心到界面的半层积分；
- $b_(i+1)$ 控制界面到第 $i+1$ 层中心的半层积分；
- 不存在由数学推导唯一给出的 $b_"eff"$。

如果把 $b_i$ 和 $b_(i+1)$ 先平均，再代入 $b / (b + 3)$，等价于假设两层之间存在一条人为构造的 Campbell 曲线。这会抹掉材料界面的水力突变，特别是在砂土/黏土相邻时会给出不可靠的界面导水率。

== 4. 更一致的界面写法

令界面水头为 $psi_I$，两侧半层厚度为：

$
  L_i = Delta z_i / 2,
  quad
  L_(i+1) = Delta z_(i+1) / 2
$

对每一侧定义 Campbell 积分原函数：

$
  F_l(psi)
  =
  integral K_l(psi) dif psi
  =
  - c_l K_l(psi) psi,
  quad
  c_l = frac(b_l, b_l + 3)
$

其中 $l$ 表示第 $i$ 层或第 $i+1$ 层。注意这里的 $c_l$ 是分层的：

$
  c_i = frac(b_i, b_i + 3),
  quad
  c_(i+1) = frac(b_(i+1), b_(i+1) + 3)
$

毛管通量连续给出：

$
  frac(F_i(psi_I) - F_i(psi_i), L_i)
  =
  frac(F_(i+1)(psi_(i+1)) - F_(i+1)(psi_I), L_(i+1))
$

这个方程求得 $psi_I$ 后，界面毛管通量为：

$
  q_psi =
  frac(F_i(psi_I) - F_i(psi_i), L_i)
  =
  frac(F_(i+1)(psi_(i+1)) - F_(i+1)(psi_I), L_(i+1))
$

再把重力项加回完整 Darcy 通量即可。这个写法的重点是：两侧积分各用各自的 $b$，界面通过通量连续连接，而不是先构造一个平均 $b$。

== 5. 不求 $psi_I$ 时的近似

若暂时不想求解界面水头 $psi_I$，也不应把 $b$ 做简单算术平均。一个比静态 $b_"eff"$ 更合理的低成本近似，是把 Campbell 修正作用在各自端点贡献上：

$
  K_(i+1/2)^"approx"
  =
  frac(
    c_i K_i psi_i
    +
    c_(i+1) K_(i+1) psi_(i+1),
    psi_i + psi_(i+1)
  )
$

其中：

$
  c_i = frac(b_i, b_i + 3),
  quad
  c_(i+1) = frac(b_(i+1), b_(i+1) + 3)
$

当 $b_i = b_(i+1) = b$ 时，上式退化为：

$
  K_(i+1/2)^"approx"
  =
  frac(b, b + 3)
  frac(K_i psi_i + K_(i+1) psi_(i+1), psi_i + psi_(i+1))
$

这与同质层的 Campbell 修正一致。异质层时，修正因子不再是一个预先平均出来的 $b_"eff"$，而是两侧端点贡献的加权结果。

若必须把它写成单个等效因子：

$
  K_(i+1/2)^"approx"
  =
  c_"eff"
  frac(K_i psi_i + K_(i+1) psi_(i+1), psi_i + psi_(i+1))
$

则：

$
  c_"eff"
  =
  frac(
    c_i K_i psi_i
    +
    c_(i+1) K_(i+1) psi_(i+1),
    K_i psi_i + K_(i+1) psi_(i+1)
  )
$

如果还要反推形式上的 $b_"eff"$，只能由：

$
  c_"eff" = frac(b_"eff", b_"eff" + 3)
$

得到：

$
  b_"eff" = frac(3 c_"eff", 1 - c_"eff")
$

这说明 $b_"eff"$ 依赖当前 $K$ 和 $psi$，不是只由 $b_i$、$b_(i+1)$ 决定的常数，更不是 $(b_i + b_(i+1)) / 2$。

== 6. 推荐结论

对于“两个层之间 $b$ 如何取值”，推荐按目标精度分三层处理：

1. 最严格：不取单一 $b$。求界面水头 $psi_I$，两侧半层分别用 $b_i$ 和 $b_(i+1)$ 做 Campbell 积分，并强制通量连续。
2. 工程近似：不求 $psi_I$，但 Campbell 修正因子分侧使用，即 $c_i = b_i / (b_i + 3)$、$c_(i+1) = b_(i+1) / (b_(i+1) + 3)$。
3. 不建议：先取 $b_"eff" = (b_i + b_(i+1)) / 2$ 再代入 $b_"eff" / (b_"eff" + 3)$。这是经验写法，不能从异质界面的 Campbell 积分推出。

因此，界面上的 $b$ 不应该被看作一个需要“平均”的参数。正确的思路是：$b$ 随材料分侧取值，界面导水率由两侧水力阻力或通量连续共同决定。
