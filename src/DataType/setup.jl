# JAX 风格的 setup 函数：将模型拆分为状态(st)和参数(ps)
#
# 设计哲学:
#   st, ps = setup(model)      # JAX/Flax 惯例
#   st = 状态变量，模拟过程中会改变
#   ps = 模型参数，模拟过程中保持不变
#
# 用法:
#   st, ps = setup(soil)                              # 从旧 Soil 转换
#   st, ps = setup(VegType, SoilType; Ta, θ0, ...)    # 直接构造
#   st, ps = setup(; VegType, SoilType, Ta, θ0, ...)  # 关键字参数版本

export setup


"""
    setup(soil::Soil; FT=Float64) -> (StateBEPS, ParamBEPS)

从传统 Soil 结构体拆分为 JAX 风格的 (state, params) 元组。

# Arguments
- `soil::Soil`: 传统的 Soil 结构体（包含状态和参数）
- `FT=Float64`: 浮点类型

# Returns
- `st::StateBEPS`: 状态变量（模拟过程中会改变）
- `ps::ParamBEPS`: 模型参数（模拟过程中保持不变）

# Example
```julia
soil = Soil()
Init_Soil_Parameters(soil, VegType, SoilType, r_root_decay)
Init_Soil_T_θ!(soil, Tsoil, Tair, θ0, z_snow)

st, ps = setup(soil)
UpdateSoilMoisture(st, ps, kstep)
```
"""
function setup(soil::Soil; FT=Float64)
  # 提取状态变量
  st = StateBEPS(soil)

  # 提取模型参数
  ps = ParamBEPS{FT}()
  Soil2Params!(ps, soil)

  return st, ps
end


"""
    setup(VegType::Int, SoilType::Int; kwargs...) -> (StateBEPS, ParamBEPS)

直接构造 JAX 风格的 (state, params) 元组，无需先创建 Soil。

# Arguments
- `VegType::Int`: 植被类型 (6=DBF, 9=EBF, 25=cropland, 40=barren, etc.)
- `SoilType::Int`: 土壤类型 (1=sand, ..., 8=silty clay loam, ..., 11=clay)

# Keyword Arguments
- `Ta::Real=20.0`: 气温 [°C]，用于初始化土壤温度
- `Tsoil::Real=Ta`: 初始土壤温度 [°C]
- `θ0::Real=0.3`: 初始土壤湿度 [m³/m³]
- `z_snow::Real=0.0`: 初始积雪深度 [m]
- `r_drainage::Real=0.5`: 地表排水速率
- `r_root_decay::Real=0.95`: 根系分布衰减率
- `N::Int=5`: 土壤层数
- `FT::Type=Float64`: 浮点类型

# Returns
- `st::StateBEPS`: 状态变量
- `ps::ParamBEPS`: 模型参数

# Example
```julia
st, ps = setup(25, 8; Ta=20.0, θ0=0.3, z_snow=0.0)
UpdateSoilMoisture(st, ps, 360.0)
```
"""
function setup(VegType::Int, SoilType::Int;
  Ta::Real=20.0,
  Tsoil::Real=Ta,
  θ0::Real=0.3,
  z_snow::Real=0.0,
  r_drainage::Real=0.5,
  r_root_decay::Real=0.95,
  N::Int=5,
  FT::Type=Float64)

  # 构造参数
  ps = ParamBEPS(VegType, SoilType; N, FT, r_drainage)
  ps.veg.r_root_decay = FT(r_root_decay)

  # 构造状态
  st = StateBEPS(; n_layer=Cint(N))

  # 复制 dz 到状态（方便计算）
  st.dz[1:N] .= ps.dz

  # 初始化根系分布
  UpdateRootFraction!(st, ps)

  # 初始化温湿度
  Init_Soil_T_θ!(st, Float64(Tsoil), Float64(Ta), Float64(θ0), Float64(z_snow))

  return st, ps
end


"""
    setup(; VegType::Int, SoilType::Int, kwargs...) -> (StateBEPS, ParamBEPS)

纯关键字参数版本，适合从 NamedTuple 或 Dict 传参。

# Example
```julia
params = (VegType=25, SoilType=8, Ta=20.0, θ0=0.3)
st, ps = setup(; params...)
```
"""
function setup(; VegType::Int, SoilType::Int, kwargs...)
  setup(VegType, SoilType; kwargs...)
end


# 泛型版本：支持传入已有的 ParamBEPS
"""
    setup(ps::ParamBEPS; kwargs...) -> (StateBEPS, ParamBEPS)

从已有参数构造状态，返回 (state, params) 元组。

# Example
```julia
ps = ParamBEPS(25, 8)
st, ps = setup(ps; Ta=20.0, θ0=0.3)
```
"""
function setup(ps::ParamBEPS;
  Ta::Real=20.0,
  Tsoil::Real=Ta,
  θ0::Real=0.3,
  z_snow::Real=0.0)

  N = ps.N
  st = StateBEPS(; n_layer=Cint(N))

  # 复制 dz 到状态
  st.dz[1:N] .= ps.dz

  # 初始化根系分布
  UpdateRootFraction!(st, ps)

  # 初始化温湿度
  Init_Soil_T_θ!(st, Float64(Tsoil), Float64(Ta), Float64(θ0), Float64(z_snow))

  return st, ps
end
