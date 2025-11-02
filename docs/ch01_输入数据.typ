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

// *预热期状态变量*： 初始化需要提供的状态变量（Optional），后续模型会自行更新迭代：

// - `Tsoil`: 可以设置为空气温度

// - `soilwater`: 70% field capacity

// - `snowdepth`: 0 mm


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

先准备驱动数据：


== 2.2 模型参数

