"""
    cal_rho_a(Tair::FT, ea::FT) where {FT<:AbstractFloat}

计算空气密度

# Arguments
- `Tair`: 气温 [°C]
- `ea`: 水汽压 [kPa]

# Returns
- 空气密度 [kg m-3]

# References
- Campbell & Norman, 1998
"""
@fastmath function cal_rho_a(Tair::FT, ea::FT) where {FT<:AbstractFloat}
  P = FT(101.3)  # [kPa] 大气压
  Rd = FT(287.0) # [J kg-1 K-1] 干空气气体常数
  T = Tair + FT(273.15)  # 转换为开尔文
  
  return P * FT(1000.0) / (Rd * T) * (1.0 - FT(0.378) * ea / P)
end

"""
    ES(T::FT) where {FT<:AbstractFloat}

计算饱和水汽压（Magnus 公式）

# Arguments
- `T`: 温度 [°C]

# Returns
- 饱和水汽压 [kPa]
"""
@fastmath function ES(T::FT) where {FT<:AbstractFloat}
  return FT(0.6108) * exp(FT(17.27) * T / (T + FT(237.3)))
end

"""
    cal_ea(Tair::FT, RH::FT) where {FT<:AbstractFloat}

计算实际水汽压

# Arguments
- `Tair`: 气温 [°C]
- `RH`: 相对湿度 [%]

# Returns
- 水汽压 [kPa]
"""
@fastmath function cal_ea(Tair::FT, RH::FT) where {FT<:AbstractFloat}
  return ES(Tair) * RH / FT(100.0)
end
