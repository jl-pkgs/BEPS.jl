function evaporation_canopy_jl(tempL::Leaf, Ta::Float64, rh_air::Float64,
  Gwater::Leaf, lai::Leaf,
  perc_water_o::Float64, perc_water_u::Float64, perc_snow_o::Float64, perc_snow_u::Float64)

  # LHw = Leaf()  # latent heat from leaves W/m2, caused by evaporation of intercepted rain
  # LHs = Leaf()  # latent heat from leaves W/m2, caused by evaporation of intercepted snow

  met = meteo_pack_jl(Ta, rh_air)
  # rho_a, cp, vpd, slope, gamma = met
  latent_water = (2.501 - 0.00237 * Ta) * 1000000
  latent_snow = 2.83 * 1000000

  # leaf level latent heat caused by evaporation or sublimation
  LHw_o_sunlit = perc_water_o * latent_heat(Ta, tempL.o_sunlit, Gwater.o_sunlit, met)
  LHw_o_shaded = perc_water_o * latent_heat(Ta, tempL.o_shaded, Gwater.o_shaded, met)  
  LHw_u_sunlit = perc_water_u * latent_heat(Ta, tempL.u_sunlit, Gwater.u_sunlit, met)
  LHw_u_shaded = perc_water_u * latent_heat(Ta, tempL.u_shaded, Gwater.u_shaded, met)

  LHs_o_sunlit = perc_snow_o * latent_heat(Ta, tempL.o_sunlit, Gwater.o_sunlit, met)
  LHs_o_shaded = perc_snow_o * latent_heat(Ta, tempL.o_shaded, Gwater.o_shaded, met)
  LHs_u_sunlit = perc_snow_u * latent_heat(Ta, tempL.u_sunlit, Gwater.u_sunlit, met)
  LHs_u_shaded = perc_snow_u * latent_heat(Ta, tempL.u_shaded, Gwater.u_shaded, met)

  evapo_water_o = 1 / (latent_water) * (LHw_o_sunlit * lai.o_sunlit + LHw_o_shaded * lai.o_shaded)
  evapo_water_u = 1 / (latent_water) * (LHw_u_sunlit * lai.u_sunlit + LHw_u_shaded * lai.u_shaded)

  evapo_snow_o = 1 / (latent_snow) * (LHs_o_sunlit * lai.o_sunlit + LHs_o_shaded * lai.o_shaded)
  evapo_snow_u = 1 / (latent_snow) * (LHs_u_sunlit * lai.u_sunlit + LHs_u_shaded * lai.u_shaded)

  return evapo_water_o, evapo_water_u, evapo_snow_o, evapo_snow_u
end

function latent_heat(Ta::Float64, Ts::Float64, gw::Float64, met::NamedTuple)
  # @unpack VPD, slope, gamma, cp, rho_a = met
  # (VPD + slope * (Ts - Ta)) * rho_a * cp * gw / gamma
  (met.VPD + met.slope * (Ts - Ta)) * met.rho_a * met.cp * gw / met.gamma
end
