
"""
    InitParam_Photo_Farquhar(VegType::Int=1; FT=Float64)

初始化 Farquhar 光合作用参数（使用硬编码默认值）

# Arguments
- `VegType::Int=1`: 植被类型代码（用于未来扩展，当前使用通用默认值）
- `FT::Type=Float64`: 浮点类型

# Returns
- `ParamPhoto_Farquhar{FT}`: 光合作用参数结构体

# Examples
```julia
params = InitParam_Photo_Farquhar(6)  # DBF（当前使用通用默认值）
```
"""
function InitParam_Photo_Farquhar(VegType::Int=1; FT=Float64)
  # 当前版本使用硬编码默认值
  # 未来可以根据 VegType 从 JSON 加载不同参数
  return ParamPhoto_Farquhar{FT}()
end


# ===== 光合作用模块参数 =====

"""
    ParamPhoto_Farquhar{FT<:AbstractFloat}

Farquhar 光合作用模型参数

# Fields
- **Vcmax25**: 25°C 时的最大羧化速率 [μmol m-2 s-1]
- **Jmax25**: 25°C 时的最大电子传递速率 [μmol m-2 s-1]
- **Rd25_ratio**: 暗呼吸与 Vcmax25 的比率 [-]
- **evc**: Vcmax 活化能 [J mol-1]
- **ejm**: Jmax 活化能 [J mol-1]
- **erd**: 暗呼吸活化能 [J mol-1]
- **toptvc**: Vcmax 最适温度 [K]
- **toptjm**: Jmax 最适温度 [K]
- **kc25**: CO2 Michaelis 常数 [μmol mol-1]
- **ko25**: O2 Michaelis 常数 [mmol mol-1]
- **tau25**: Rubisco 特异性因子 [mmol mol-1]
- **qalpha**: 量子效率 [-]
- **theta2**: 电子传输曲率参数 [-]
- **g0_w**: Ball-Berry 气孔导度截距 [mol m-2 s-1]
- **g1_w**: Ball-Berry 气孔导度斜率 [-]
- **clumping**: 叶片聚集指数 [-]
"""
@bounds @with_kw mutable struct ParamPhoto_Farquhar{FT<:AbstractFloat}
  # ===== Farquhar 模型核心参数 =====
  Vcmax25::FT = 89.45 | (5.0, 200.0)              # [μmol m-2 s-1] 最大羧化速率
  Jmax25::FT = 2.39 * 89.45 - 14.2 | (50.0, 400.0) # [μmol m-2 s-1] 最大电子传递速率

  # ===== 呼吸参数 =====
  Rd25_ratio::FT = 0.004657 | (0.001, 0.01)       # [-] Rd/Vcmax25 比率
  Rd_light_factor::FT = 0.4 | (0.2, 0.6)          # [-] 光下呼吸降低因子

  # ===== 温度响应参数 [J mol-1] =====
  evc::FT = 30000.0 | (20000.0, 80000.0)          # Vcmax 活化能
  ejm::FT = 40000.0 | (20000.0, 80000.0)          # Jmax 活化能
  erd::FT = 38000.0 | (20000.0, 60000.0)          # 暗呼吸活化能

  # ===== 最适温度 [K] =====
  toptvc::FT = 298.15 | (288.15, 308.15)          # Vcmax 最适温度
  toptjm::FT = 298.15 | (288.15, 308.15)          # Jmax 最适温度

  # ===== Michaelis-Menten 常数 =====
  kc25::FT = 404.0 | (200.0, 600.0)               # [μmol mol-1] CO2
  ko25::FT = 248000.0 | (200000.0, 300000.0)      # [mmol mol-1] O2
  tau25::FT = 2600.0 | (2000.0, 3000.0)           # [mmol mol-1] Rubisco 特异性

  # ===== 电子传输参数 =====
  qalpha::FT = 0.3 | (0.2, 0.5)                   # [-] 量子效率
  theta2::FT = 0.7 | (0.5, 0.9)                   # [-] 曲率参数

  # ===== 气孔导度参数 =====
  g0_w::FT = 0.001 | (0.0, 0.01)                  # [mol m-2 s-1] Ball-Berry 截距
  g1_w::FT = 4.0 | (1.0, 10.0)                    # [-] Ball-Berry 斜率

  # ===== 冠层结构参数 =====
  clumping::FT = 0.85 | (0.5, 1.0)                # [-] 叶片聚集指数
end

export ParamPhoto_Farquhar
export InitParam_Photo_Farquhar
