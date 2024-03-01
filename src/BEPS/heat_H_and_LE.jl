function latent_heat(Ta::Float64, Ts::Float64, gw::Float64, met::NamedTuple)
  (; VPD, Δ, γ, cp, ρₐ) = met
  (VPD + Δ * (Ts - Ta)) * ρₐ * cp * gw / γ
  # (met.VPD + met.slope * (Ts - Ta)) * met.ρₐ * met.cp * gw / met.γ
end

function latent_heat!(leleaf::Leaf, Gw::Leaf, VPD, slope, Tc_old::Leaf, Tair, ρₐ, cp, γ)
  leleaf.o_sunlit = Gw.o_sunlit * (VPD + slope * (Tc_old.o_sunlit - Tair)) * ρₐ * cp / γ
  leleaf.o_shaded = Gw.o_shaded * (VPD + slope * (Tc_old.o_shaded - Tair)) * ρₐ * cp / γ
  leleaf.u_sunlit = Gw.u_sunlit * (VPD + slope * (Tc_old.u_sunlit - Tair)) * ρₐ * cp / γ
  leleaf.u_shaded = Gw.u_shaded * (VPD + slope * (Tc_old.u_shaded - Tair)) * ρₐ * cp / γ
end

function transpiration_jl(T_leaf::Leaf, Ta::Float64, RH::Float64, Gtrans::Leaf, lai::Leaf)
  T = Leaf() # transpiration

  met = meteo_pack_jl(Ta, RH)
  (; ρₐ, cp, VPD, Δ, γ) = met
  λ = cal_lambda(Ta)

  # Luo, 2018, JGR-Biogeosciences
  T.o_sunlit = (VPD + Δ * (T_leaf.o_sunlit - Ta)) * ρₐ * cp * Gtrans.o_sunlit / γ
  T.o_shaded = (VPD + Δ * (T_leaf.o_shaded - Ta)) * ρₐ * cp * Gtrans.o_shaded / γ
  T.u_sunlit = (VPD + Δ * (T_leaf.u_sunlit - Ta)) * ρₐ * cp * Gtrans.u_sunlit / γ
  T.u_shaded = (VPD + Δ * (T_leaf.u_shaded - Ta)) * ρₐ * cp * Gtrans.u_shaded / γ

  trans_o = 1.0 / λ * (T.o_sunlit * lai.o_sunlit + T.o_shaded * lai.o_shaded)
  trans_u = 1.0 / λ * (T.u_sunlit * lai.u_sunlit + T.u_shaded * lai.u_shaded)

  trans_o, trans_u
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
  (; ρₐ, cp) = met # cp: specific heat of moist air above canopy
  
  SH = sensible_heat(T_leaf, Ta, ρₐ, cp, Gheat)
  SH_o::FT = SH.o_sunlit * lai.o_sunlit + SH.o_shaded * lai.o_shaded
  SH_u::FT = SH.u_sunlit * lai.u_sunlit + SH.u_shaded * lai.u_shaded

  SH_o = max(-200.0, SH_o)
  SH_u = max(-200.0, SH_u)
  SH_g::FT = (T_ground - Ta) * ρₐ * cp * Gheat_g

  SH_o, SH_u, SH_g
end
