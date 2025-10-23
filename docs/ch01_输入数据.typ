#import "@local/modern-cug-report:0.1.3": *
#show: doc => template(doc, footer: "BESP", header: "1. 输入数据")


= 1 *模型使用*

模型运行需要输入气象驱动（forcing），叶面积指数LAI，以及部分土壤和植被参数。

```julia
forcing = deserialize("data/p1_meteo")
lai = readdlm("examples/input/p1_lai.txt")[:]

par = (lon=114.3197, lat=30.35, clumping=0.75,
  landcover=13,
  soil_type=7,
  Tsoil=4.2, soilwater=0.203, snowdepth=0.0)
```

*预热期状态变量*： 初始化需要提供的状态变量（Optional），后续模型会自行更新迭代：

- `Tsoil`: 可以设置为空气温度

- `soilwater`: 70% field capacity

- `snowdepth`: 0 mm



= 2 *输入数据*

== 2.1 气象驱动

- `doy`: day of year

- `hour`: hour of day, UTC时间

- `rad`: 入射短波辐射，Rs_in, W m-2

- `tem`: 气温，Tair, ℃

- `hum`: 比湿，q (specific humidity), g/kg

- `pre`: 降水，prcp (precipitation), mm

- `wind`: 2m风速，U2, m/s

*不需要输入长波辐射数据。*

== 2.2 土壤参数

// *BEPS模型中使用的土壤类型编号：*

// #figure(
//   caption: [],
//   table(
//     columns: 5,
//     // (2.0cm, 1.9cm, 1.55cm),
//     rows: (0.8cm, 0.8cm),
//     align: (horizon, horizon + left, horizon, horizon, horizon+left),
    
//     [*Soil ID*], [*Soil Name*], [], [*Soil ID*], [*Soil Name*],
//     [1], [sand], [], [2], [loamy sand],
//     [3], [sandy loam], [], [4], [loam],
//     [5], [silt loam], [], [6], [sandy clay loam],
//     [7], [clay loam], [], [8], [silty clay loam],
//     [9], [sandy clay], [], [10], [silty clay],
//     [11], [clay], [], [12 (default)], [silt],
//   ),
// ) <table_>

#figure(
  caption: [BEPS土壤类型。],
  table(
    columns: 6,
    // (2.0cm, 1.9cm, 1.55cm),
    rows: (0.8cm, 0.77cm),
    align: (horizon),

    [*Value*], [*Color Value*], [*Description*], [*Longname*], [*BEPS*], [*Name*],
    [1], [\#d5c36b], [Cl], [clay], [11], [粘土],
    [2], [\#b96947], [SiCl], [silty clay], [10], [粉砂质粘土],
    [3], [\#9d3706], [SaCl], [sandy clay], [9], [砂质粘土],
    [4], [\#ae868f], [ClLo], [clay loam], [7], [粘土壤],
    [5], [\#f86714], [SiClLo], [silty clay loam], [8], [粉砂质粘土壤],
    [6], [\#46d143], [SaClLo], [sandy clay loam], [6], [砂质粘土壤],
    [7], [\#368f20], [Lo], [loam], [4], [壤土],
    [8], [\#3e5a14], [SiLo], [silt loam], [5], [粉壤土],
    [9], [\#ffd557], [SaLo], [sandy loam], [3], [砂壤土],
    [10], [\#fff72e], [Si], [silt], [12 (default)], [粉壤土],
    [11], [\#ff5a9d], [LoSa], [loamy sand], [2], [粘土壤],
    [12], [\#ff005b], [Sa], [sand], [1], [砂土],
  ),
) <table_SoilType>

#table-note[
  *注*：表格第一列Value，采用的是OpenLandMap Soil Texture Class (USDA System)土壤类型编号。
]

提供土壤类型编号，是最简单的土壤参数输入方式，模型会根据土壤类型编号，自动调用相应的土壤参数。
*如果精细化地控制各个参数，可参考后续章节提供中间变量。*

土壤驱动数据，可从GEE下载，https://code.earthengine.google.com/2b14ec8d599d7dcaa1dedac47b261aca?noload=1。

=== 2.2.1 Optional variables:

- `d_soil`       : 每层土壤厚度，the depth of each soil layer, (m)
- `density_soil` : 土壤密度（干土），soil bulk density, (kg/m3)
- `f_org`        : 土壤孔隙度，volume fraction of organic matter

  #beamer-block[用于模拟土壤温度与土壤湿度。]

=== 2.2.2 土壤水力与热传导参数

- `b`            : Cambell parameter b
- `Ksat`         : 饱和水力传导系数，saturated hydraulic conductivity
- `fei`          : 孔隙度，porosity
- `theta_vfc`    : 田间持水量，field capacity (ignored)
- `θ_vwp`        : 凋萎含水量，wilt point
- `thermal_cond` : 热传导系数$kappa$，thermal conductivity
- `psi_sat`      : 饱和土水势$psi_"sat"$，water potential at saturated


== 2.3 植被类型参数

- 冠层聚集度$Omega$：He et al., (2012)

#figure(
  caption: [BEPS植被类型参数],
  table(
    columns: 7,
    // (2.0cm, 1.9cm, 1.55cm),
    rows: (0.8cm, 0.8cm),
    align: (horizon),
    [_*IGBP*_],
    [*$italic("LC")$*],
    [*$italic("LC")_{"code"}$*],
    [],
    [_*IGBP*_],
    [*$italic("LC")$*],
    [*$italic("LC")_{"code"}$*],

    [DBF], [DBF], [6], [], [EBF], [EBF], [9],
    [ENF], [ENF], [1], [], [DNF], [DNF], [2],
    [CRO], [Others], [-1], [], [CSH], [Shrub], [13],
    [OSH], [Shrub], [13], [], [SAV], [Shrub], [13],
    [WSA], [Shrub], [13], [], [MF], [Shrub], [13],
    [GRA], [Others], [-1], [], [SNO], [Others], [-1],
    [WET], [Others], [-1],
  ),
) <table_VegType>

#table-note[
  注：其中$italic("LC")$、$italic("LC")_{"code"}$为BEPS模型中采用的植被类型与植被类型编号。`IGBP`可从MOD12Q1数据集中获取。
]
