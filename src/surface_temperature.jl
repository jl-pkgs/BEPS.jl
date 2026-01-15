# 计算雪热导率 (Jordan, 1991)
@inline cal_κ_snow(ρ) = 0.021 + 4.2e-4 * ρ + 2.2e-9 * ρ^3

# 显式时间步进 (Explicit Euler)
@inline step_exp(T, Fin, Fout, C, dt) = T + (Fin - Fout) / C * dt

# 相变检查：若跨越0度且存在水源，则强制钳制在0度
@inline check_phase(T, Told, w) = ((T > 0 >= Told) || (T < 0 <= Told && w)) ? zero(T) : T

# Formula: (T_old*Inertia + Rad + Air + Cond_Below) / Denom
function solve_imp(T_old, T_bnd, T_bot, ΔE, ra, z, G, ρCp, κ_bot; z_rad=z, c_s=1.0, μ=NaN)
  I = ΔE * ra * z
  numerator = T_old * I + G * ra * z_rad + ρCp * T_bnd * z + c_s * ra * κ_bot * T_bot
  denominator = ρCp * z + c_s * ra * κ_bot + I

  ans = numerator / denominator
  !isnan(μ) && (ans = clamp(ans, μ - 25.0, μ + 25.0))
  return ans
end


function surface_temperature_jl!(
  Rn_g::FT, T_air::FT, Tc_u::FT, RH::FT, 
  z_snow::FT, z_water::FT, ρ_snow::FT, perc_snow_g::FT,
  z_soil1::FT, κ_soil1::FT, Cv_soil1::FT, Cv_soil0::FT, Gheat_g::FT,
  E_soil::FT, E_water_g::FT, E_snow_g::FT,
  G_soil1::FT, T_soil1_last::FT, T_soil0_last::FT,
  last::SnowLand{FT}, current::SnowLand{FT}
) where {FT<:AbstractFloat}
  # 从 last_in 读取上一步的温度
  T_surf_last = last.T_surf
  T_mix0_last = last.T_mix0
  T_snow0_last = last.T_snow0
  T_snow1_last = last.T_snow1
  T_snow2_last = last.T_snow2

  Δt::FT = kstep
  cp_ice::FT = 2228.261
  λ_water::FT = cal_lambda(T_air)
  λ_snow::FT = 2.83e6
  cp::FT = cal_cp(T_air, RH)
  ra_g::FT = 1.0 / Gheat_g
  ρCp::FT = ρₐ * cp

  # Parameters for Case 3 (Deep Snow)
  dz_snow_s1::FT = 0.02
  dz_snow_s2::FT = 0.02
  dz_soil_s0::FT = 0.02

  κ_dry_snow::FT = cal_κ_snow(ρ_snow) # 雪热导率
  Gg::FT = Rn_g - E_snow_g * λ_snow - (E_water_g + E_soil) * λ_water # 地表可用能量

  T_soil0::FT = 0.0
  G::FT = 0.0

  ΔM_soil1 = Cv_soil1 * 0.02 / Δt # soil heat capacity per unit area, Cv = ρ cp

  if z_snow <= 0.02
    # Case 1: 无雪或极浅雪 (≤2cm)
    T_surf = solve_imp(T_surf_last, T_air, T_soil1_last, ΔM_soil1,
      ra_g, z_soil1, Gg, ρCp, κ_soil1; μ=T_surf_last)

    T_mix0 = T_surf
    T_snow0 = T_mix0
    T_soil0 = T_mix0
    T_snow1 = T_mix0
    T_snow2 = T_mix0

    G = 2 * κ_soil1 * (T_mix0 - T_soil1_last) / z_soil1
    G = clamp(G, -100.0, 100.0)

  elseif z_snow > 0.02 && z_snow <= 0.05
    # Case 2: 中等雪深 (2-5cm) - 雪土混合
    Δz_soil1 = 0.5z_soil1 # 第一层土壤厚度的一半
    Δz_snow = z_snow      # 雪层厚度, bottom to top

    # Case 2: 中等雪深 (2-5cm) - 雪土混合
    T_soil0 = solve_imp(T_soil0_last, T_air, T_soil1_last, ΔM_soil1, ra_g, z_soil1,
      Gg, ρCp, κ_soil1; c_s=2.0, μ=T_air)                                       # 裸土地表温度

    ΔM_snow = cp_ice * ρ_snow * z_snow / Δt
    T_snow0 = solve_imp(T_snow0_last, Tc_u, T_mix0_last, ΔM_snow, ra_g, z_snow,
      Gg, ρCp, κ_dry_snow; μ=T_air)                                             # 雪表温度

    T_interface = (κ_soil1 * T_soil1_last / Δz_soil1 + κ_dry_snow * T_snow0 / Δz_snow + ΔM_soil1 * T_mix0_last) /
                  (κ_soil1 / Δz_soil1 + κ_dry_snow / Δz_snow + ΔM_soil1)        # Ts(z=0) only with snow cover
    T_mix0 = T_soil0 * (1 - perc_snow_g) + T_interface * perc_snow_g       # Ts(z=0) 

    G_snow = κ_dry_snow * (T_snow0 - T_soil1_last) / (Δz_snow + Δz_soil1)
    G_soil = κ_soil1 * (T_mix0 - T_soil1_last) / Δz_soil1

    G = G_snow * perc_snow_g + G_soil * (1 - perc_snow_g)
    G = clamp(G, -100.0, 100.0)

    # 相变检查
    (T_snow0 > 0.0 && T_snow0_last <= 0.0 && z_snow > 0.0) && (T_snow0 = 0.0)
    (T_snow0 < 0.0 && T_snow0_last >= 0.0 && z_water > 0.0) && (T_snow0 = 0.0)

    T_surf = T_snow0 * perc_snow_g + T_soil0 * (1 - perc_snow_g)
    T_surf = clamp(T_surf, T_air - 25.0, T_air + 25.0)

    T_snow1 = T_snow0
    T_snow2 = T_snow0

  else  # z_snow > 0.05
    # Case 3: 深雪 (>5cm) - 3层雪模型
    dz_snow_s12 = dz_snow_s1 + dz_snow_s2

    ΔM_snow = cp_ice * ρ_snow * dz_snow_s1 / Δt
    T_snow0 = solve_imp(T_snow0_last, T_air, T_snow1_last, ΔM_snow, ra_g, dz_snow_s12, Gg, ρCp, κ_dry_snow; z_rad=dz_snow_s1, μ=T_air)

    G_snow = κ_dry_snow * (T_snow0 - T_snow1_last) / dz_snow_s12
    G = clamp(G_snow, -100.0, 100.0)

    G_snow1 = κ_dry_snow * (T_snow1_last - T_snow2_last) / (z_snow - dz_snow_s1)
    T_snow1 = step_exp(T_snow1_last, G, G_snow1, cp_ice * ρ_snow * dz_snow_s2, Δt)

    G_snow2 = (T_snow2_last - T_mix0_last) / ((0.5 * (z_snow - dz_snow_s12) / κ_dry_snow) + (dz_soil_s0 / κ_soil1))
    T_snow2 = step_exp(T_snow2_last, G_snow1, G_snow2, cp_ice * ρ_snow * (z_snow - dz_snow_s12), Δt)

    T_mix0 = step_exp(T_mix0_last, G_snow2, G_soil1, Cv_soil0 * dz_soil_s0, Δt)
    T_soil0 = T_mix0

    # 相变检查
    (T_snow0 > 0.0 && T_snow0_last <= 0.0 && z_snow > 0.0) && (T_snow0 = 0.0)
    (T_snow0 < 0.0 && T_snow0_last >= 0.0 && z_water > 0.0) && (T_snow0 = 0.0)

    T_surf = T_snow0
  end
  @pack! current = T_surf, T_mix0, T_snow0, T_snow1, T_snow2
  return G, T_soil0
end


function surface_temperature!(
  soil::AbstractSoil, Tsnow_p::SnowLand{FT}, Tsnow_c::SnowLand{FT}, 
  radiation_g::FT, Tc_u::FT, T_air::FT, RH::FT, z_snow::FT, z_water::FT, ρ_snow::FT, f_snow_g::FT, Gheat_g::FT,
  Evap_soil::FT, Evap_SW::FT, Evap_SS::FT) where {FT<:AbstractFloat}

  # 1. 准备土壤热力学参数
  Cv1 = Cv0 = soil.Cv[1]           # 体积热容 [J m-3 K-1]
  κ_soil1 = soil.κ[1]                 # 热导率 [W m-1 K-1]; 表皮为1, 第一层为2
  z_soil1 = soil.dz[1]               # 层厚度 [m]

  # 2. 调用 surface_temperature_jl 计算
  G, T_soil0 = surface_temperature_jl!(radiation_g, T_air, Tc_u, RH,
    z_snow, z_water, ρ_snow, f_snow_g, 
    z_soil1, κ_soil1, Cv1, Cv0, Gheat_g,
    Evap_soil, Evap_SW, Evap_SS,
    soil.G[1],
    soil.Tsoil_p[2],      # T_soil1: 第一层土壤温度 [°C]
    soil.Tsoil_p[1],      # T_soil0:
    Tsnow_p, Tsnow_c)

  soil.Tsoil_c[1] = T_soil0 # 更新 soil 的当前温度状态
  return G
end
