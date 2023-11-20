# 1. 土壤参数

土壤驱动数据，可从[GEE](https://code.earthengine.google.com/2b14ec8d599d7dcaa1dedac47b261aca?noload=1)下载。

## 1.1. 方案1：提供土壤类型(`soil_type`)

- 1: `sand`
- 2: `loamy sand`
- 3: `sandy loam`
- 4: `loam`
- 5: `silt loam`
- 6: `sandy clay loam`
- 7: `clay loam`
- 8: `silty clay loam`
- 9: `sandy clay`
- 10: `silty clay`
- 11: `clay`
- default:


## 1.2. 方案2：提供中间变量数据

### 1.2.1. Optional variables:

- `d_soil`       : the depth of each soil layer, (m)
- `density_soil` : soil bulk density, (kg/m3)
- `f_org`        : volume fraction of organic matter

### 1.2.2. Hard variables: 
- `b`            : Cambell parameter b
- `Ksat`         : saturated hydraulic conductivity
- `fei`          : porosity
- `theta_vfc`    : field capacity (ignored)
- `theta_vwp`    : wilt point
- `thermal_cond` : thermal conductivity
- `psi_sat`      : water potential at saturated


# 2. 其他土壤参数

- `Tsoil`: 可以设置为空气温度

- `soilwater`: 70% field capacity

- `snowdepth`: 0 mm
