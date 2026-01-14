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


function surface_temperature_jl(T_air::FT, RH::FT, z_snow::FT, z_water::FT,
  Cv_soil1::FT, Cv_soil0::FT, Gheat_g::FT,
  z_soil1::FT, ρ_snow::FT, tempL_u::FT, Rn_g::FT,
  E_soil::FT, E_water_g::FT, E_snow_g::FT, κ_soil1::FT,
  perc_snow_g::FT, G_soil1::FT,
  T_ground_last::FT,
  T_soil1_last::FT, T_soil_surf_last::FT, T_soil0_last::FT,
  T_snow_last::FT, T_snow1_last::FT, T_snow2_last::FT) where {FT<:AbstractFloat}

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

  T_ground::FT = 0.0
  T_soil_surf::FT = 0.0
  T_soil0::FT = 0.0
  T_snow::FT = 0.0
  T_snow1::FT = 0.0
  T_snow2::FT = 0.0
  G::FT = 0.0

  ΔM_soil1 = Cv_soil1 * 0.02 / Δt # soil heat capacity per unit area, Cv = ρ cp

  if z_snow <= 0.02
    # Case 1: 无雪或极浅雪 (≤2cm)
    T_ground = solve_imp(T_ground_last, T_air, T_soil1_last, ΔM_soil1, 
      ra_g, z_soil1, Gg, ρCp, κ_soil1; μ=T_ground_last)

    T_soil_surf = T_ground
    T_snow = T_soil_surf
    T_soil0 = T_soil_surf
    T_snow1 = T_soil_surf
    T_snow2 = T_soil_surf

    G = 2 * κ_soil1 * (T_soil_surf - T_soil1_last) / z_soil1
    G = clamp(G, -100.0, 100.0)

  elseif z_snow > 0.02 && z_snow <= 0.05
    # Case 2: 中等雪深 (2-5cm) - 雪土混合
    Δz_soil1 = 0.5z_soil1 # 第一层土壤厚度的一半
    Δz_snow = z_snow      # 雪层厚度, bottom to top

    # Case 2: 中等雪深 (2-5cm) - 雪土混合
    T_soil0 = solve_imp(T_soil0_last, T_air, T_soil1_last, ΔM_soil1, ra_g, z_soil1, 
      Gg, ρCp, κ_soil1; c_s=2.0, μ=T_air)                                       # 裸土地表温度

    ΔM_snow = cp_ice * ρ_snow * z_snow / Δt
    T_snow = solve_imp(T_snow_last, tempL_u, T_soil_surf_last, ΔM_snow, ra_g, z_snow, 
      Gg, ρCp, κ_dry_snow; μ=T_air)                                             # 雪表温度

    T_interface = (κ_soil1 * T_soil1_last / Δz_soil1 + κ_dry_snow * T_snow / Δz_snow + ΔM_soil1 * T_soil_surf_last) /
                  (κ_soil1 / Δz_soil1 + κ_dry_snow / Δz_snow + ΔM_soil1)        # Ts(z=0) only with snow cover
    T_soil_surf = T_soil0 * (1 - perc_snow_g) + T_interface * perc_snow_g       # Ts(z=0) 

    G_snow = κ_dry_snow * (T_snow - T_soil1_last) / (Δz_snow + Δz_soil1)
    G_soil = κ_soil1 * (T_soil_surf - T_soil1_last) / Δz_soil1

    G = G_snow * perc_snow_g + G_soil * (1 - perc_snow_g)
    G = clamp(G, -100.0, 100.0)

    # 相变检查
    (T_snow > 0.0 && T_snow_last <= 0.0 && z_snow > 0.0) && (T_snow = 0.0)
    (T_snow < 0.0 && T_snow_last >= 0.0 && z_water > 0.0) && (T_snow = 0.0)

    T_ground = T_snow * perc_snow_g + T_soil0 * (1 - perc_snow_g)
    T_ground = clamp(T_ground, T_air - 25.0, T_air + 25.0)

    T_snow1 = T_snow
    T_snow2 = T_snow

  else  # z_snow > 0.05
    # Case 3: 深雪 (>5cm) - 3层雪模型
    dz_snow_s12 = dz_snow_s1 + dz_snow_s2
    
    ΔM_snow = cp_ice * ρ_snow * dz_snow_s1 / Δt
    T_snow = solve_imp(T_snow_last, T_air, T_snow1_last, ΔM_snow, ra_g, dz_snow_s12, Gg, ρCp, κ_dry_snow; z_rad=dz_snow_s1, μ=T_air)

    G_snow = κ_dry_snow * (T_snow - T_snow1_last) / dz_snow_s12
    G = clamp(G_snow, -100.0, 100.0)

    G_snow1 = κ_dry_snow * (T_snow1_last - T_snow2_last) / (z_snow - dz_snow_s1)
    T_snow1 = step_exp(T_snow1_last, G, G_snow1, cp_ice * ρ_snow * dz_snow_s2, Δt)

    G_snow2 = (T_snow2_last - T_soil_surf_last) / ((0.5 * (z_snow - dz_snow_s12) / κ_dry_snow) + (dz_soil_s0 / κ_soil1))
    T_snow2 = step_exp(T_snow2_last, G_snow1, G_snow2, cp_ice * ρ_snow * (z_snow - dz_snow_s12), Δt)

    T_soil_surf = step_exp(T_soil_surf_last, G_snow2, G_soil1, Cv_soil0 * dz_soil_s0, Δt)
    T_soil0 = T_soil_surf

    # 相变检查
    (T_snow > 0.0 && T_snow_last <= 0.0 && z_snow > 0.0) && (T_snow = 0.0)
    (T_snow < 0.0 && T_snow_last >= 0.0 && z_water > 0.0) && (T_snow = 0.0)

    T_ground = T_snow
  end
  return G, T_ground, T_soil_surf, T_soil0, T_snow, T_snow1, T_snow2
end



"""
    surface_temperature!(soil, cache, k, T_air, RH,
        z_snow, z_water, Gheat_g, ρ_snow, Tc_u, radiation_g,
        Evap_soil, Evap_SW, Evap_SS, f_snow_g, dz, κ)

辅助函数：准备数据并计算表面温度

# 功能说明
这个函数封装了以下操作：
1. 从 soil 提取并准备热力学参数（Cs, κ, dz, G）
2. 调用 surface_temperature_jl 计算表面温度
3. 将结果保存到 cache（用于子时间步循环）和 soil

# 参数
- `soil`: 土壤状态结构体
- `cache`: 中间临时变量结构体
- `k`: 当前时间步索引
- 其他参数同 surface_temperature_jl

# 副作用
- 修改 cache.G[1], cache.T_ground[k], cache.T_surf_mix[k],
  cache.T_surf_snow[k], cache.T_snow_L1[k], cache.T_snow_L2[k]
- 修改 soil.Tsoil_c[1]
- 修改 dz[2], κ[2]（临时数组）
"""
function surface_temperature!(
  soil::AbstractSoil, cache::TransientCache, k::Int,
  T_air::FT, RH::FT, z_snow::FT, z_water::FT,
  Gheat_g::FT, ρ_snow::FT, Tc_u::FT,
  radiation_g::FT, Evap_soil::FT, Evap_SW::FT, Evap_SS::FT,
  f_snow_g::FT, dz::Vector{FT}, κ::Vector{FT}) where {FT<:AbstractFloat}

  # 1. 准备土壤热力学参数
  cache.Cs[1:2] .= soil.Cs[1]      # 体积热容 [J m-3 K-1]
  cache.Tc_u[k] = Tc_u                # 下层冠层温度 [°C]
  κ[2] = soil.κ[1]                    # 热导率 [W m-1 K-1]; 表皮为1, 第一层为2
  dz[2] = soil.dz[1]                  # 层厚度 [m]
  cache.G[2] = soil.G[1]           # 第1层土壤热通量 [W m-2]

  # 2. 调用 surface_temperature_jl 计算
  # 直接使用 soil.Tsoil_p 作为输入，避免通过 cache.T_soil 中转
  G, T_ground, T_soil_surf, T_soil0, T_snow, T_snow1, T_snow2 =
    surface_temperature_jl(T_air, RH, z_snow, z_water,
      cache.Cs[2], cache.Cs[1], Gheat_g, dz[2], ρ_snow, cache.Tc_u[k],
      radiation_g, Evap_soil, Evap_SW, Evap_SS,
      κ[2], f_snow_g, cache.G[2],
      cache.T_ground[k-1],
      soil.Tsoil_p[2],      # T_soil1_last: 第一层土壤温度 [°C]
      soil.Tsoil_p[1],      # T_soil_surf_last: 土壤表面温度（雪-土交界面）[°C]
      cache.T_surf_mix[k-1],
      cache.T_surf_snow[k-1], cache.T_snow_L1[k-1], cache.T_snow_L2[k-1])

  # 3. 保存结果到 cache（用于子时间步循环）
  cache.G[1] = G
  cache.T_ground[k] = T_ground
  cache.T_surf_mix[k] = T_soil0
  cache.T_surf_snow[k] = T_snow
  cache.T_snow_L1[k] = T_snow1
  cache.T_snow_L2[k] = T_snow2

  # 4. 更新 soil 的当前温度状态
  soil.Tsoil_c[1] = T_soil_surf
  return nothing
end
