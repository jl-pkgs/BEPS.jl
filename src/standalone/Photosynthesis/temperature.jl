"""
    TBOLTZ(T::FT, E::FT) where {FT<:AbstractFloat}

温度响应函数（Arrhenius 方程）

# Arguments
- `T`: 温度 [K]
- `E`: 活化能 [J mol-1]

# Returns
- 温度响应因子 [0-1]

# References
- Campbell & Norman, 1998
"""
@fastmath function TBOLTZ(T::FT, E::FT) where {FT<:AbstractFloat}
  rugc = FT(8.314)  # [J mol-1 K-1] 通用气体常数
  TK25 = FT(298.16) # [K] 25°C
  return exp((T - TK25) * E / (T * TK25 * rugc))
end

"""
    fTv(T::FT, Vcmax25::FT, evc::FT, toptvc::FT) where {FT<:AbstractFloat}

Vcmax 的温度响应（Medlyn et al. 2002 归一化峰值 Arrhenius）

# Arguments
- `T`: 叶片温度 [K]
- `Vcmax25`: 25°C 时的最大羧化速率 [μmol m-2 s-1]
- `evc`: 活化能 [J mol-1]
- `toptvc`: 最适温度 [K]（当前实现未使用，保留以保持 API 兼容性；
  去活化由固定参数 Hd=200000 J/mol, S=640 J/mol/K 控制）

# Returns
- 温度调整后的 Vcmax [μmol m-2 s-1]，在 25°C 归一化为 Vcmax25

# References
- Medlyn et al., 2002, Plant Cell Environ.
"""
@fastmath function fTv(T::FT, Vcmax25::FT, evc::FT, toptvc::FT) where {FT<:AbstractFloat}
  rugc = FT(8.314)
  TK25 = FT(298.15)
  Hd   = FT(200000.0)  # 去活化焓 [J mol-1]
  S    = FT(640.0)     # 熵项 [J mol-1 K-1]

  f_temp  = TBOLTZ(T, evc)
  num = FT(1.0) + exp((S * TK25 - Hd) / (TK25 * rugc))  # 25°C 归一化因子
  den = FT(1.0) + exp((S * T - Hd) / (T * rugc))          # 当前温度去活化
  return Vcmax25 * f_temp * num / den
end

"""
    fTj(T::FT, Jmax25::FT, ejm::FT, toptjm::FT) where {FT<:AbstractFloat}

Jmax 的温度响应（Medlyn et al. 2002 归一化峰值 Arrhenius）

# Arguments
- `T`: 叶片温度 [K]
- `Jmax25`: 25°C 时的最大电子传递速率 [μmol m-2 s-1]
- `ejm`: 活化能 [J mol-1]
- `toptjm`: 最适温度 [K]（当前实现未使用，保留以保持 API 兼容性；
  去活化由固定参数 Hd=200000 J/mol, S=640 J/mol/K 控制）

# Returns
- 温度调整后的 Jmax [μmol m-2 s-1]，在 25°C 归一化为 Jmax25

# References
- Medlyn et al., 2002, Plant Cell Environ.
"""
@fastmath function fTj(T::FT, Jmax25::FT, ejm::FT, toptjm::FT) where {FT<:AbstractFloat}
  rugc = FT(8.314)
  TK25 = FT(298.15)
  Hd   = FT(200000.0)  # 去活化焓 [J mol-1]
  S    = FT(640.0)     # 熵项 [J mol-1 K-1]

  f_temp  = TBOLTZ(T, ejm)
  num = FT(1.0) + exp((S * TK25 - Hd) / (TK25 * rugc))  # 25°C 归一化因子
  den = FT(1.0) + exp((S * T - Hd) / (T * rugc))          # 当前温度去活化
  return Jmax25 * f_temp * num / den
end

"""
    fTd(T::FT, Rd25::FT, erd::FT) where {FT<:AbstractFloat}

暗呼吸的温度响应

# Arguments
- `T`: 叶片温度 [K]
- `Rd25`: 25°C 时的暗呼吸速率 [μmol m-2 s-1]
- `erd`: 活化能 [J mol-1]

# Returns
- 温度调整后的 Rd [μmol m-2 s-1]
"""
@fastmath function fTd(T::FT, Rd25::FT, erd::FT) where {FT<:AbstractFloat}
  return Rd25 * TBOLTZ(T, erd)
end
