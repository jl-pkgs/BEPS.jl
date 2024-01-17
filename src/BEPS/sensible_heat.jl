function sensible_heat_jl(T_leaf::Leaf, T_ground::FT, Ta::FT, RH::FT,
  Gheat::Leaf, Gheat_g::FT, lai::Leaf)

  SH = Leaf()
  met = meteo_pack_jl(Ta, RH)
  ρₐ = met.rho_a
  cp = met.cp     # specific heat of moist air above canopy

  SH.o_sunlit = (T_leaf.o_sunlit - Ta) * ρₐ * cp * Gheat.o_sunlit
  SH.o_shaded = (T_leaf.o_shaded - Ta) * ρₐ * cp * Gheat.o_shaded
  SH.u_sunlit = (T_leaf.u_sunlit - Ta) * ρₐ * cp * Gheat.u_sunlit
  SH.u_shaded = (T_leaf.u_shaded - Ta) * ρₐ * cp * Gheat.u_shaded

  SH_o::FT = SH.o_sunlit * lai.o_sunlit + SH.o_shaded * lai.o_shaded
  SH_u::FT = SH.u_sunlit * lai.u_sunlit + SH.u_shaded * lai.u_shaded

  SH_o = max(-200.0, SH_o)
  SH_u = max(-200.0, SH_u)
  SH_g::FT = (T_ground - Ta) * ρₐ * cp * Gheat_g

  SH_o, SH_u, SH_g
end


function sensible_heat(T_w::FT, T_a::FT, ρₐ::FT, cp::FT, gH::FT)::FT
  (T_w - T_a) * ρₐ * cp * gH
end

function sensible_heat(T_leaf::Leaf, T_a::FT, ρₐ::FT, cp::FT, gH::Leaf)::FT
  SH = Leaf()
  
  SH.o_sunlit = (T_leaf.o_sunlit - T_a) * ρₐ * cp * gH.o_sunlit
  SH.o_shaded = (T_leaf.o_shaded - T_a) * ρₐ * cp * gH.o_shaded
  SH.u_sunlit = (T_leaf.u_sunlit - T_a) * ρₐ * cp * gH.u_sunlit
  SH.u_shaded = (T_leaf.u_shaded - T_a) * ρₐ * cp * gH.u_shaded
  # SH.o_sunlit = sensible_heat(T_leaf.o_sunlit, gH.o_sunlit, T_a, ρₐ, cp)
  # SH.o_shaded = sensible_heat(T_leaf.o_shaded, gH.o_shaded, T_a, ρₐ, cp)
  # SH.u_sunlit = sensible_heat(T_leaf.u_sunlit, gH.u_sunlit, T_a, ρₐ, cp)
  # SH.u_shaded = sensible_heat(T_leaf.u_shaded, gH.u_shaded, T_a, ρₐ, cp)
  SH
end
