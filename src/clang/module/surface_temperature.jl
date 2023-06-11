# function surface_temperature(T_air, rh_air, depth_snow, depth_water, capacity_heat_soil1, capacity_heat_soil0, Gheat_g, depth_soil1, density_snow, tempL_u, netRad_g, evapo_soil, evapo_water_g, evapo_snow_g, lambda_soil1, percent_snow_g, heat_flux_soil1, 
#   T_ground_last, T_soil1_last, T_any0_last, T_snow_last, T_soil0_last, T_snow1_last, T_snow2_last,
#   T_ground::TypeRef, T_any0::TypeRef, T_snow::TypeRef, T_soil0::TypeRef, 
#   T_snow1::TypeRef, T_snow2::TypeRef, heat_flux::TypeRef)

#   ccall((:surface_temperature, libbeps), Cvoid,
#     (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble,
#       Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
#     T_air, rh_air, depth_snow, depth_water, capacity_heat_soil1, capacity_heat_soil0, Gheat_g, depth_soil1, density_snow, tempL_u, netRad_g, evapo_soil, evapo_water_g, evapo_snow_g, lambda_soil1, percent_snow_g, heat_flux_soil1, T_ground_last, T_soil1_last, T_any0_last, T_snow_last, T_soil0_last, T_snow1_last, T_snow2_last,
#     T_ground, T_any0, T_snow, T_soil0, T_snow1, T_snow2, heat_flux)
# end
function surface_temperature(T_air, rh_air, depth_snow, depth_water, capacity_heat_soil1, capacity_heat_soil0, Gheat_g, depth_soil1, density_snow, tempL_u, netRad_g, evapo_soil, evapo_water_g, evapo_snow_g, lambda_soil1, percent_snow_g, heat_flux_soil1, T_ground_last, T_soil1_last, T_any0_last, T_snow_last, T_soil0_last, T_snow1_last, T_snow2_last)

  T_ground = Ref(0.0)
  T_any0 = Ref(0.0)
  T_snow = Ref(0.0)
  T_soil0 = Ref(0.0)
  T_snow1 = Ref(0.0)
  T_snow2 = Ref(0.0)
  heat_flux = Ref(0.0)

  ccall((:surface_temperature, libbeps), Cvoid,
    (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble,
      Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
    T_air, rh_air, depth_snow, depth_water, capacity_heat_soil1, capacity_heat_soil0, Gheat_g, depth_soil1, density_snow, tempL_u, netRad_g, evapo_soil, evapo_water_g, evapo_snow_g, lambda_soil1, percent_snow_g, heat_flux_soil1, T_ground_last, T_soil1_last, T_any0_last, T_snow_last, T_soil0_last, T_snow1_last, T_snow2_last, T_ground, T_any0, T_snow, T_soil0, T_snow1, T_snow2, heat_flux)

  T_ground[], T_any0[], T_snow[], T_soil0[], T_snow1[], T_snow2[], heat_flux[]
end


function surface_temperature_jl(T_air::T, rh_air::T, depth_snow::T, depth_water::T,
  capacity_heat_soil1::T, capacity_heat_soil0::T, Gheat_g::T,
  depth_soil1::T, density_snow::T, tempL_u::T, netRad_g::T,
  evapo_soil::T, evapo_water_g::T, evapo_snow_g::T, lambda_soil1::T,
  percent_snow_g::T, heat_flux_soil1::T, T_ground_last::T,
  T_soil1_last::T, T_any0_last::T, T_snow_last::T,
  T_soil0_last::T, T_snow1_last::T, T_snow2_last::T) where {T<:Real}

  length_step = kstep
  cp_ice = 2228.261  # specific heat of ice
  latent_water = cal_lambda(T_air) # J kg-1
  latent_snow = 2.83 * 1000000

  density_air = rho_a
  cp_air = cal_cp(T_air, rh_air)

  ra_g = 1 / Gheat_g  # aerodynamic resistance of heat

  # thermal conductivity of snow
  lambda_snow = 0.021 + 4.2 * density_snow / 10000 + 2.2 * density_snow^3 * 1e-9

  # available energy on ground for
  Gg = netRad_g - evapo_snow_g * latent_snow - (evapo_water_g + evapo_soil) * latent_water

  T_ground, T_any0, T_snow, T_soil0, T_snow1, T_snow2 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
  heat_flux = 0.0

  if depth_snow <= 0.02
    ttt = capacity_heat_soil1 * 0.02 / length_step
    T_ground = (T_ground_last * ttt * ra_g * depth_soil1 + Gg * ra_g * depth_soil1 + density_air * cp_air * T_air * depth_soil1 + ra_g * lambda_soil1 * T_soil1_last) /
                  (density_air * cp_air * depth_soil1 + ra_g * lambda_soil1 + ttt * ra_g * depth_soil1)
    T_ground = clamp(T_ground, T_ground_last - 25, T_ground_last + 25)

    T_any0 = T_ground
    T_snow = T_any0
    T_soil0 = T_any0
    T_snow1 = T_any0
    T_snow2 = T_any0

    heat_flux = 2 * lambda_soil1 * (T_any0 - T_soil1_last) / depth_soil1
    heat_flux = clamp(heat_flux, -100, 100)
  elseif depth_snow > 0.02 && depth_snow <= 0.05
    ttt = capacity_heat_soil1 * 0.02 / length_step  # for soil fraction part

    T_soil0 = (T_soil0_last * ttt * ra_g * depth_soil1 + Gg * ra_g * depth_soil1 + density_air * cp_air * T_air * depth_soil1 + 2 * ra_g * lambda_soil1 * T_soil1_last) /
                 (density_air * cp_air * depth_soil1 + 2 * ra_g * lambda_soil1 + ttt * ra_g * depth_soil1)

    T_soil0 = clamp(T_soil0, T_air - 25, T_air + 25)

    ttt = cp_ice * density_snow * depth_snow / length_step  # for snow part
    T_snow = (T_snow_last * ttt * ra_g * depth_snow + Gg * ra_g * depth_snow + density_air * cp_air * tempL_u * depth_snow + ra_g * lambda_snow * T_any0_last) /
                (density_air * cp_air * depth_snow + ra_g * lambda_snow + ttt * ra_g * depth_snow)

    T_snow = clamp(T_snow, T_air - 25, T_air + 25)

    ttt = (lambda_soil1 * T_soil1_last / depth_soil1 + T_snow * lambda_snow + 0.02 * capacity_heat_soil1 / length_step * T_any0_last) /
          (lambda_soil1 / depth_soil1 + lambda_snow / depth_snow + 0.02 * capacity_heat_soil1 / length_step)
    T_any0 = T_soil0 * (1 - percent_snow_g) + ttt * percent_snow_g

    heat_flux_snow = lambda_snow / (depth_snow + 0.5 * depth_soil1) * (T_snow - T_soil1_last)
    heat_flux_soil = heat_flux_snow * (T_any0 - T_soil1_last) / depth_soil1

    heat_flux = heat_flux_snow * percent_snow_g + heat_flux_soil * (1 - percent_snow_g)
    heat_flux = clamp(heat_flux, -100, 100)

    # starting to melt
    if T_snow > zero && T_snow_last <= zero && depth_snow > zero
      T_snow = 0
    end

    # starting to frozen
    if T_snow < zero && T_snow_last >= zero && depth_water > zero
      T_snow = 0
    end

    # percent_snow_g =min(1.0,Wg_snow[kkk] / (0.05 * rho_snow[kkk])); # use the fraction before
    T_ground = T_snow * percent_snow_g + T_soil0 * (1 - percent_snow_g)
    T_ground = clamp(T_ground, T_air - 25, T_air + 25)

    T_snow1 = T_snow
    T_snow2 = T_snow
  elseif depth_snow > 0.05
    ttt = cp_ice * density_snow * 0.02 / length_step

    T_snow = (T_snow_last * ttt * ra_g * 0.04 + Gg * ra_g * 0.02 + density_air * cp_air * T_air * 0.04 + ra_g * lambda_snow * T_snow1_last) /
                (density_air * cp_air * 0.04 + ra_g * lambda_snow + ttt * ra_g * 0.04)
    T_snow = clamp(T_snow, T_air - 25, T_air + 25)

    heat_flux_snow = lambda_snow * (T_snow - T_snow1_last) / 0.04

    heat_flux = heat_flux_snow
    heat_flux = clamp(heat_flux, -100, 100)

    heat_flux_snow1 = lambda_snow * (T_snow1_last - T_snow2_last) / (depth_snow - 0.02)
    T_snow1 = T_snow1_last + ((heat_flux - heat_flux_snow1) / (cp_ice * density_snow * 0.02)) * length_step
    heat_flux_snow2 = (T_snow2_last - T_any0_last) / ((0.5 * (depth_snow - 0.04) / lambda_snow) + (0.02 / lambda_soil1))
    T_snow2 = T_snow2_last + ((heat_flux_snow1 - heat_flux_snow2) / (cp_ice * density_snow * (depth_snow - 0.04))) * length_step

    T_any0 = T_any0_last + ((heat_flux_snow2 - heat_flux_soil1) / (capacity_heat_soil0 * 0.02)) * length_step
    T_soil0 = T_any0

    if T_snow > zero && T_snow_last <= zero && depth_snow > zero
      T_snow = 0
    end

    if T_snow < zero && T_snow_last >= zero && depth_water > zero
      T_snow = 0
    end
    T_ground = T_snow
  end

  return T_ground, T_any0, T_snow, T_soil0, T_snow1, T_snow2, heat_flux
end
