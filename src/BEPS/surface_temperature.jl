function surface_temperature_jl(T_air::FT, rh_air::FT, depth_snow::FT, depth_water::FT,
  capacity_heat_soil1::FT, capacity_heat_soil0::FT, Gheat_g::FT,
  depth_soil1::FT, ρ_snow::FT, tempL_u::FT, Rn_g::FT,
  evapo_soil::FT, evapo_water_g::FT, evapo_snow_g::FT, λ_soil1::FT,
  perc_snow_g::FT, G_soil1::FT, 
  T_ground_last::FT,
  T_soil1_last::FT, T_any0_last::FT, T_soil0_last::FT,
  T_snow_last::FT, T_snow1_last::FT, T_snow2_last::FT)

  length_step::FT = kstep
  cp_ice::FT = 2228.261  # specific heat of ice
  latent_water::FT = cal_lambda(T_air) # J kg-1
  latent_snow::FT = 2.83 * 1000000

  ρ_air::FT = rho_a
  cp::FT = cal_cp(T_air, rh_air)

  ra_g::FT = 1.0 / Gheat_g  # aerodynamic resistance of heat

  # thermal conductivity of snow
  λ_snow::FT = 0.021 + 4.2 * ρ_snow / 10000 + 2.2 * ρ_snow^3 * 1e-9

  # available energy on ground for
  Gg::FT = Rn_g - evapo_snow_g * latent_snow - (evapo_water_g + evapo_soil) * latent_water

  T_ground::FT = 0.0
  T_any0::FT = 0.0
  T_soil0::FT = 0.0
  T_snow::FT = 0.0
  T_snow1::FT = 0.0
  T_snow2::FT = 0.0
  G::FT = 0.0

  if depth_snow <= 0.02
    ttt = capacity_heat_soil1 * 0.02 / length_step
    T_ground = (T_ground_last * ttt * ra_g * depth_soil1 + Gg * ra_g * depth_soil1 + ρ_air * cp * T_air * depth_soil1 + ra_g * λ_soil1 * T_soil1_last) /
               (ρ_air * cp * depth_soil1 + ra_g * λ_soil1 + ttt * ra_g * depth_soil1)
    T_ground = clamp(T_ground, T_ground_last - 25, T_ground_last + 25)

    T_any0 = T_ground
    T_snow = T_any0
    T_soil0 = T_any0
    T_snow1 = T_any0
    T_snow2 = T_any0

    G = 2 * λ_soil1 * (T_any0 - T_soil1_last) / depth_soil1
    G = clamp(G, -100.0, 100.0)

  elseif depth_snow > 0.02 && depth_snow <= 0.05
    ttt = capacity_heat_soil1 * 0.02 / length_step  # for soil fraction part

    T_soil0 = (T_soil0_last * ttt * ra_g * depth_soil1 + Gg * ra_g * depth_soil1 + ρ_air * cp * T_air * depth_soil1 + 2 * ra_g * λ_soil1 * T_soil1_last) /
              (ρ_air * cp * depth_soil1 + 2 * ra_g * λ_soil1 + ttt * ra_g * depth_soil1)

    T_soil0 = clamp(T_soil0, T_air - 25.0, T_air + 25.0)

    ttt = cp_ice * ρ_snow * depth_snow / length_step  # for snow part
    T_snow = (T_snow_last * ttt * ra_g * depth_snow + Gg * ra_g * depth_snow + ρ_air * cp * tempL_u * depth_snow + ra_g * λ_snow * T_any0_last) /
             (ρ_air * cp * depth_snow + ra_g * λ_snow + ttt * ra_g * depth_snow)

    T_snow = clamp(T_snow, T_air - 25.0, T_air + 25.0)

    ttt = (λ_soil1 * T_soil1_last / depth_soil1 + T_snow * λ_snow + 0.02 * capacity_heat_soil1 / length_step * T_any0_last) /
          (λ_soil1 / depth_soil1 + λ_snow / depth_snow + 0.02 * capacity_heat_soil1 / length_step)
    T_any0 = T_soil0 * (1 - perc_snow_g) + ttt * perc_snow_g

    G_snow = λ_snow / (depth_snow + 0.5 * depth_soil1) * (T_snow - T_soil1_last)
    G_soil = G_snow * (T_any0 - T_soil1_last) / depth_soil1

    G = G_snow * perc_snow_g + G_soil * (1 - perc_snow_g)
    G = clamp(G, -100.0, 100.0)

    # starting to melt
    if T_snow > zero && T_snow_last <= zero && depth_snow > zero
      T_snow = 0.0
    end

    # starting to frozen
    if T_snow < zero && T_snow_last >= zero && depth_water > zero
      T_snow = 0.0
    end

    # perc_snow_g =min(1.0,Wg_snow[kkk] / (0.05 * rho_snow[kkk])); # use the fraction before
    T_ground = T_snow * perc_snow_g + T_soil0 * (1 - perc_snow_g)
    T_ground = clamp(T_ground, T_air - 25.0, T_air + 25.0)

    T_snow1 = T_snow
    T_snow2 = T_snow
  elseif depth_snow > 0.05
    ttt = cp_ice * ρ_snow * 0.02 / length_step

    T_snow = (T_snow_last * ttt * ra_g * 0.04 + Gg * ra_g * 0.02 + ρ_air * cp * T_air * 0.04 + ra_g * λ_snow * T_snow1_last) /
             (ρ_air * cp * 0.04 + ra_g * λ_snow + ttt * ra_g * 0.04)
    T_snow = clamp(T_snow, T_air - 25.0, T_air + 25.0)

    G_snow = λ_snow * (T_snow - T_snow1_last) / 0.04
    G = clamp(G_snow, -100.0, 100.0)

    G_snow1 = λ_snow * (T_snow1_last - T_snow2_last) / (depth_snow - 0.02)
    T_snow1 = T_snow1_last + ((G - G_snow1) / (cp_ice * ρ_snow * 0.02)) * length_step
    G_snow2 = (T_snow2_last - T_any0_last) / ((0.5 * (depth_snow - 0.04) / λ_snow) + (0.02 / λ_soil1))
    T_snow2 = T_snow2_last + ((G_snow1 - G_snow2) / (cp_ice * ρ_snow * (depth_snow - 0.04))) * length_step

    T_any0 = T_any0_last + ((G_snow2 - G_soil1) / (capacity_heat_soil0 * 0.02)) * length_step
    T_soil0 = T_any0

    if T_snow > zero && T_snow_last <= zero && depth_snow > zero
      T_snow = 0
    end

    if T_snow < zero && T_snow_last >= zero && depth_water > zero
      T_snow = 0
    end
    T_ground = T_snow
  end

  return G, T_ground, T_any0, T_soil0, T_snow, T_snow1, T_snow2
end

