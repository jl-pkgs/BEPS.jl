"""
    surface_temperature_jl(T_air, rh_air, z_snow, z_water,
      cp_soil1, cp_soil0, Gheat_g,
      z_soil1, ρ_snow, tempL_u, Rn_g,
      E_soil, E_water_g, E_snow_g, λ_soil1,
      perc_snow_g, G_soil1,
      T_ground_last,
      T_soil1_last, T_soil_surf_last, T_soil0_last,
      T_snow_last, T_snow1_last, T_snow2_last)

计算地表温度和土壤热通量。

# 状态变量 (上一时间步) (State Variables)
- `T_ground_last`: t-1 时刻的混合地表温度 [°C]。指积雪表面（如果有雪）或土壤表面的温度。
- `T_soil1_last`: t-1 时刻的第一层土壤温度 [°C]。
- `T_soil_surf_last`: t-1 时刻的土壤表面温度 [°C]。指雪-土交界面温度（如果有雪）或土壤表面温度。
- `T_soil0_last`: t-1 时刻的表层土壤等效温度 [°C]。深雪模型的中间变量。
- `T_snow_last`: t-1 时刻的积雪表面温度 [°C]。
- `T_snow1_last`: t-1 时刻的积雪第一层（上层）温度 [°C]（用于深雪模型）。
- `T_snow2_last`: t-1 时刻的积雪第二层（下层）温度 [°C]（用于深雪模型）。

# 返回值 (Returns)
- `G`: 地表热通量 [W m-2]
- `T_ground`: 混合地表温度 [°C]
- `T_soil_surf`: 土壤表面温度 [°C]
- `T_soil0`: 表层土壤等效温度 [°C]
- `T_snow`: 积雪表面温度 [°C]
- `T_snow1`: 积雪第一层温度 [°C]
- `T_snow2`: 积雪第二层温度 [°C]

# 示意图 (Diagram)
```text
       Air (T_air)
          |
          | ra_g (Aerodynamic Resistance)
          |
    ------------------  <-- T_ground (Bulk Surface Temperature)
    |   Snow Layer   |      (If snow exists, T_ground ≈ T_snow)
    |  [ T_snow1 ]   |  <-- T_snow1 (Upper Snow Layer Temp)
    |----------------|
    |  [ T_snow2 ]   |  <-- T_snow2 (Lower Snow Layer Temp)
    ------------------  <-- T_soil_surf (Soil Surface Temp / Snow-Soil Interface)
          |                 (If no snow, T_soil_surf ≈ T_ground)
          |
      Soil Layer 0      <-- T_soil0 (Effective Surface Soil Temp)
          |
    ------------------
          |
      Soil Layer 1      <-- T_soil1 (Layer 1 Soil Temp)
          |
```
"""
function surface_temperature_jl(T_air::FT, rh_air::FT, z_snow::FT, z_water::FT,
  cp_soil1::FT, cp_soil0::FT, Gheat_g::FT,
  z_soil1::FT, ρ_snow::FT, tempL_u::FT, Rn_g::FT,
  E_soil::FT, E_water_g::FT, E_snow_g::FT, λ_soil1::FT,
  perc_snow_g::FT, G_soil1::FT,
  T_ground_last::FT,
  T_soil1_last::FT, T_soil_surf_last::FT, T_soil0_last::FT,
  T_snow_last::FT, T_snow1_last::FT, T_snow2_last::FT)

  Δt::FT = kstep
  cp_ice::FT = 2228.261  # specific heat of ice
  λ_water::FT = cal_lambda(T_air) # J kg-1
  λ_snow::FT = 2.83 * 1000000

  cp::FT = cal_cp(T_air, rh_air)
  ra_g::FT = 1.0 / Gheat_g  # aerodynamic resistance of heat

  # thermal conductivity of snow
  κ_dry_snow::FT = 0.021 + 4.2 * ρ_snow / 10000 + 2.2 * ρ_snow^3 * 1e-9

  # available energy on ground for
  Gg::FT = Rn_g - E_snow_g * λ_snow - (E_water_g + E_soil) * λ_water

  T_ground::FT = 0.0
  T_soil_surf::FT = 0.0
  T_soil0::FT = 0.0
  T_snow::FT = 0.0
  T_snow1::FT = 0.0
  T_snow2::FT = 0.0
  G::FT = 0.0

  if z_snow <= 0.02
    ΔE = cp_soil1 * 0.02 / Δt
    T_ground = (T_ground_last * ΔE * ra_g * z_soil1 + Gg * ra_g * z_soil1 + ρₐ * cp * T_air * z_soil1 + ra_g * λ_soil1 * T_soil1_last) /
               (ρₐ * cp * z_soil1 + ra_g * λ_soil1 + ΔE * ra_g * z_soil1)
    T_ground = clamp(T_ground, T_ground_last - 25, T_ground_last + 25)

    T_soil_surf = T_ground
    T_snow = T_soil_surf
    T_soil0 = T_soil_surf
    T_snow1 = T_soil_surf
    T_snow2 = T_soil_surf

    G = 2 * λ_soil1 * (T_soil_surf - T_soil1_last) / z_soil1
    G = clamp(G, -100.0, 100.0)

  elseif z_snow > 0.02 && z_snow <= 0.05
    ΔE_soil = cp_soil1 * 0.02 / Δt  # for soil fraction part

    T_soil0 = (T_soil0_last * ΔE_soil * ra_g * z_soil1 + Gg * ra_g * z_soil1 + ρₐ * cp * T_air * z_soil1 + 2 * ra_g * λ_soil1 * T_soil1_last) /
              (ρₐ * cp * z_soil1 + 2 * ra_g * λ_soil1 + ΔE_soil * ra_g * z_soil1)

    T_soil0 = clamp(T_soil0, T_air - 25.0, T_air + 25.0)

    ΔE_snow = cp_ice * ρ_snow * z_snow / Δt  # for snow fraction part
    T_snow = (T_snow_last * ΔE_snow * ra_g * z_snow + Gg * ra_g * z_snow + ρₐ * cp * tempL_u * z_snow + ra_g * κ_dry_snow * T_soil_surf_last) /
             (ρₐ * cp * z_snow + ra_g * κ_dry_snow + ΔE_snow * ra_g * z_snow)

    T_snow = clamp(T_snow, T_air - 25.0, T_air + 25.0)

    T_weighted = (λ_soil1 * T_soil1_last / z_soil1 + T_snow * κ_dry_snow + 0.02 * cp_soil1 / Δt * T_soil_surf_last) /
                 (λ_soil1 / z_soil1 + κ_dry_snow / z_snow + 0.02 * cp_soil1 / Δt)
    T_soil_surf = T_soil0 * (1 - perc_snow_g) + T_weighted * perc_snow_g

    G_snow = κ_dry_snow / (z_snow + 0.5 * z_soil1) * (T_snow - T_soil1_last)
    G_soil = G_snow * (T_soil_surf - T_soil1_last) / z_soil1

    G = G_snow * perc_snow_g + G_soil * (1 - perc_snow_g)
    G = clamp(G, -100.0, 100.0)

    # starting to melt
    if T_snow > 0.0 && T_snow_last <= 0.0 && z_snow > 0.0
      T_snow = 0.0
    end
    # starting to frozen
    if T_snow < 0.0 && T_snow_last >= 0.0 && z_water > 0.0
      T_snow = 0.0
    end

    # perc_snow_g =min(1.0,Wcs_g[kkk] / (0.05 * rho_snow[kkk])); # use the fraction before
    T_ground = T_snow * perc_snow_g + T_soil0 * (1 - perc_snow_g)
    T_ground = clamp(T_ground, T_air - 25.0, T_air + 25.0)

    T_snow1 = T_snow
    T_snow2 = T_snow
  elseif z_snow > 0.05
    ΔE_snow = cp_ice * ρ_snow * 0.02 / Δt

    T_snow = (T_snow_last * ΔE_snow * ra_g * 0.04 + Gg * ra_g * 0.02 + ρₐ * cp * T_air * 0.04 + ra_g * κ_dry_snow * T_snow1_last) /
             (ρₐ * cp * 0.04 + ra_g * κ_dry_snow + ΔE_snow * ra_g * 0.04)
    T_snow = clamp(T_snow, T_air - 25.0, T_air + 25.0)

    G_snow = κ_dry_snow * (T_snow - T_snow1_last) / 0.04
    G = clamp(G_snow, -100.0, 100.0)

    G_snow1 = κ_dry_snow * (T_snow1_last - T_snow2_last) / (z_snow - 0.02)
    T_snow1 = T_snow1_last + ((G - G_snow1) / (cp_ice * ρ_snow * 0.02)) * Δt
    G_snow2 = (T_snow2_last - T_soil_surf_last) / ((0.5 * (z_snow - 0.04) / κ_dry_snow) + (0.02 / λ_soil1))
    T_snow2 = T_snow2_last + ((G_snow1 - G_snow2) / (cp_ice * ρ_snow * (z_snow - 0.04))) * Δt

    T_soil_surf = T_soil_surf_last + ((G_snow2 - G_soil1) / (cp_soil0 * 0.02)) * Δt
    T_soil0 = T_soil_surf

    if T_snow > 0.0 && T_snow_last <= 0.0 && z_snow > 0.0
      T_snow = 0
    end
    if T_snow < 0.0 && T_snow_last >= 0.0 && z_water > 0.0
      T_snow = 0
    end
    T_ground = T_snow
  end

  return G, T_ground, T_soil_surf, T_soil0, T_snow, T_snow1, T_snow2
end

