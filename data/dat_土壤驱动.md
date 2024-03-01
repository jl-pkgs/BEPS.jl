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
- default (12): silt

**表1. USDA土壤类型分类**

| Value | Color Value | Description | Longname        | BEPS |
| :---- | :---------- | :---------- | --------------- | ---- |
| 1     | #d5c36b     | Cl          | clay            | 11   |
| 2     | #b96947     | SiCl        | silty clay      | 10   |
| 3     | #9d3706     | SaCl        | sandy clay      | 9    |
| 4     | #ae868f     | ClLo        | clay loam       | 7    |
| 5     | #f86714     | SiClLo      | silty clay loam | 8    |
| 6     | #46d143     | SaClLo      | sandy clay loam | 6    |
| 7     | #368f20     | Lo          | loam            | 4    |
| 8     | #3e5a14     | SiLo        | silt loam       | 5    |
| 9     | #ffd557     | SaLo        | sandy loam      | 3    |
| 10    | #fff72e     | Si          | silt            | 12   |
| 11    | #ff5a9d     | LoSa        | loamy sand      | 2    |
| 12    | #ff005b     | Sa          | sand            | 1    |


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
- `θ_vwp`        : wilt point
- `thermal_cond` : thermal conductivity
- `psi_sat`      : water potential at saturated


# 2. 其他土壤参数

- `Tsoil`: 可以设置为空气温度

- `soilwater`: 70% field capacity

- `snowdepth`: 0 mm
