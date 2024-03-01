function latent_heat(Ta::Float64, Ts::Float64, gw::Float64, met::NamedTuple)
  # @unpack VPD, slope, gamma, cp, ρ_a = met
  # (VPD + slope * (Ts - Ta)) * ρ_a * cp * gw / gamma
  (met.VPD + met.slope * (Ts - Ta)) * met.ρ_a * met.cp * gw / met.gamma
end

function latent_heat!(leleaf::Leaf, Gw::Leaf, VPD, slope, Tc_old::Leaf, Tair, ρ_a, Cp_ca, gamma)
  leleaf.o_sunlit = Gw.o_sunlit * (VPD + slope * (Tc_old.o_sunlit - Tair)) * ρ_a * Cp_ca / gamma
  leleaf.o_shaded = Gw.o_shaded * (VPD + slope * (Tc_old.o_shaded - Tair)) * ρ_a * Cp_ca / gamma
  leleaf.u_sunlit = Gw.u_sunlit * (VPD + slope * (Tc_old.u_sunlit - Tair)) * ρ_a * Cp_ca / gamma
  leleaf.u_shaded = Gw.u_shaded * (VPD + slope * (Tc_old.u_shaded - Tair)) * ρ_a * Cp_ca / gamma
end


function sensible_heat(T_w::FT, T_a::FT, ρₐ::FT, cp::FT, gH::FT)::FT
  (T_w - T_a) * ρₐ * cp * gH
end

function sensible_heat(T_leaf::Leaf, T_a::FT, ρₐ::FT, cp::FT, gH::Leaf)
  SH = Leaf()  
  SH.o_sunlit = (T_leaf.o_sunlit - T_a) * ρₐ * cp * gH.o_sunlit
  SH.o_shaded = (T_leaf.o_shaded - T_a) * ρₐ * cp * gH.o_shaded
  SH.u_sunlit = (T_leaf.u_sunlit - T_a) * ρₐ * cp * gH.u_sunlit
  SH.u_shaded = (T_leaf.u_shaded - T_a) * ρₐ * cp * gH.u_shaded

  SH
end

function sensible_heat_jl(T_leaf::Leaf, T_ground::FT, Ta::FT, RH::FT,
  Gheat::Leaf, Gheat_g::FT, lai::Leaf)

  met = meteo_pack_jl(Ta, RH)
  ρₐ = met.ρ_a
  cp = met.cp     # specific heat of moist air above canopy

  SH = sensible_heat(T_leaf, Ta, ρₐ, cp, Gheat)
  
  SH_o::FT = SH.o_sunlit * lai.o_sunlit + SH.o_shaded * lai.o_shaded
  SH_u::FT = SH.u_sunlit * lai.u_sunlit + SH.u_shaded * lai.u_shaded

  SH_o = max(-200.0, SH_o)
  SH_u = max(-200.0, SH_u)
  SH_g::FT = (T_ground - Ta) * ρₐ * cp * Gheat_g

  SH_o, SH_u, SH_g
end
