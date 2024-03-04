function evaporation_canopy_jl(T_leaf::Leaf, Ta::Float64, RH::Float64,
  Gwater::Leaf, lai::Leaf,
  perc_water::Layer2{Float64}, 
  perc_snow::Layer3{Float64})
  # perc_water_o::Float64, perc_water_u::Float64, 
  # perc_snow_o::Float64, perc_snow_u::Float64)

  # LHw = Leaf()  # latent heat from leaves W/m2, caused by evaporation of intercepted rain
  # LHs = Leaf()  # latent heat from leaves W/m2, caused by evaporation of intercepted snow
  met = meteo_pack_jl(Ta, RH)
  λ = met.λ # 2.5*e6 J / kg
  
  # leaf level latent heat caused by evaporation or sublimation
  LHw_o_sunlit = perc_water.o * latent_heat(Ta, T_leaf.o_sunlit, Gwater.o_sunlit, met)
  LHw_o_shaded = perc_water.o * latent_heat(Ta, T_leaf.o_shaded, Gwater.o_shaded, met)  
  LHw_u_sunlit = perc_water.u * latent_heat(Ta, T_leaf.u_sunlit, Gwater.u_sunlit, met)
  LHw_u_shaded = perc_water.u * latent_heat(Ta, T_leaf.u_shaded, Gwater.u_shaded, met)

  LHs_o_sunlit = perc_snow.o * latent_heat(Ta, T_leaf.o_sunlit, Gwater.o_sunlit, met)
  LHs_o_shaded = perc_snow.o * latent_heat(Ta, T_leaf.o_shaded, Gwater.o_shaded, met)
  LHs_u_sunlit = perc_snow.u * latent_heat(Ta, T_leaf.u_sunlit, Gwater.u_sunlit, met)
  LHs_u_shaded = perc_snow.u * latent_heat(Ta, T_leaf.u_shaded, Gwater.u_shaded, met)

  E_water_o = (LHw_o_sunlit * lai.o_sunlit + LHw_o_shaded * lai.o_shaded) / λ
  E_water_u = (LHw_u_sunlit * lai.u_sunlit + LHw_u_shaded * lai.u_shaded) / λ

  E_snow_o = (LHs_o_sunlit * lai.o_sunlit + LHs_o_shaded * lai.o_shaded) / λ_snow
  E_snow_u = (LHs_u_sunlit * lai.u_sunlit + LHs_u_shaded * lai.u_shaded) / λ_snow

  return E_water_o, E_water_u, E_snow_o, E_snow_u
end
