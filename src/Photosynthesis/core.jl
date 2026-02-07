using UnPack

"""
    farquhar_model(T::FT, PAR::FT, ci::FT, params::ParamPhoto_Farquhar{FT}) where {FT<:AbstractFloat}

单叶片 Farquhar 光合作用模型

# Arguments
- `T`: 叶片温度 [K]
- `PAR`: 光合有效辐射 [μmol m-2 s-1]
- `ci`: 胞间 CO2 浓度 [μmol mol-1]
- `params`: 光合作用参数

# Returns
- `An`: 净光合速率 [μmol m-2 s-1]
- `Rd`: 暗呼吸速率 [μmol m-2 s-1]

# References
- Farquhar et al., 1980
"""
function farquhar_model(T::FT, PAR::FT, ci::FT, params) where {FT<:AbstractFloat}
  @unpack Vcmax25, Jmax25, Rd25_ratio, Rd_light_factor,
  evc, ejm, erd, toptvc, toptjm,
  kc25, ko25, tau25, qalpha, theta2 = params

  # 1. 温度调整
  Vcmax = fTv(T, Vcmax25, evc, toptvc)
  Jmax = fTj(T, Jmax25, ejm, toptjm)
  Rd25 = Rd25_ratio * Vcmax25
  Rd = fTd(T, Rd25, erd)

  # 2. Michaelis-Menten 常数
  Kc = kc25 * TBOLTZ(T, FT(80500.0))
  Ko = ko25 * TBOLTZ(T, FT(14500.0))
  tau = tau25 * TBOLTZ(T, FT(-29000.0))

  # 3. Rubisco 限制速率
  # Wc = Vcmax * (ci - Γ*) / (ci + Kc * (1 + O2/Ko))
  Γstar = 40.0 / tau  # CO2 补偿点
  Wc = Vcmax * (ci - Γstar) / (ci + Kc * (1.0 + 210.0 / Ko))

  # 4. 光限制速率
  # J = (αI + Jmax - sqrt((αI + Jmax)^2 - 4θαIJmax)) / 2θ
  alpha = qalpha
  I2 = alpha * PAR
  J = (I2 + Jmax - sqrt((I2 + Jmax)^2 - 4.0 * theta2 * I2 * Jmax)) / (2.0 * theta2)
  Wj = J * (ci - Γstar) / (4.0 * ci + 8.0 * Γstar)

  # 5. TPU 限制（简化）
  Wp = FT(0.167) * Vcmax

  # 6. 净光合速率
  Ac = min(Wc, Wj, Wp)
  An = Ac - Rd

  return An, Rd
end
