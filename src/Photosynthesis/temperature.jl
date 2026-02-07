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

Vcmax 的温度响应

# Arguments
- `T`: 叶片温度 [K]
- `Vcmax25`: 25°C 时的最大羧化速率 [μmol m-2 s-1]
- `evc`: 活化能 [J mol-1]
- `toptvc`: 最适温度 [K]

# Returns
- 温度调整后的 Vcmax [μmol m-2 s-1]
"""
@fastmath function fTv(T::FT, Vcmax25::FT, evc::FT, toptvc::FT) where {FT<:AbstractFloat}
  rugc = FT(8.314)
  
  # 温度响应因子
  ft = TBOLTZ(T, evc)
  
  # 高温抑制因子
  fth = 1.0 + exp((-22000.0 + 710.0 * (T - 298.15)) / (298.15 * rugc))
  fth2 = 1.0 + exp((200000.0 - 640.0 * T) / (T * rugc))
  
  return Vcmax25 * ft / fth * fth2
end

"""
    fTj(T::FT, Jmax25::FT, ejm::FT, toptjm::FT) where {FT<:AbstractFloat}

Jmax 的温度响应

# Arguments
- `T`: 叶片温度 [K]
- `Jmax25`: 25°C 时的最大电子传递速率 [μmol m-2 s-1]
- `ejm`: 活化能 [J mol-1]
- `toptjm`: 最适温度 [K]

# Returns
- 温度调整后的 Jmax [μmol m-2 s-1]
"""
@fastmath function fTj(T::FT, Jmax25::FT, ejm::FT, toptjm::FT) where {FT<:AbstractFloat}
  rugc = FT(8.314)
  
  ft = TBOLTZ(T, ejm)
  fth = 1.0 + exp((30000.0 - 630.0 * (T - 298.15)) / (298.15 * rugc))
  fth2 = 1.0 + exp((100000.0 - 200.0 * T) / (T * rugc))
  
  return Jmax25 * ft / fth * fth2
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
