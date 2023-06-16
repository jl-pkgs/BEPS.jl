function sensible_heat_c(tempL::Leaf, 
  temp_g::Cdouble, temp_air::Cdouble, RH::Cdouble,
  Gheat::Leaf, Gheat_g::Cdouble, LAI::Leaf)

  SH_o = init_dbl()
  SH_u = init_dbl()
  SH_g = init_dbl()

  ccall((:sensible_heat, libbeps), Cvoid,
    (Leaf, Cdouble, Cdouble, Cdouble, Leaf, Cdouble, Leaf,
      Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
    tempL, temp_g, temp_air, RH, Gheat, Gheat_g, LAI,
    SH_o, SH_u, SH_g)
  
  SH_o[], SH_u[], SH_g[]
end


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


function sensible_heat(T_w::FT, T_a::FT, gH::FT, ρ_a::FT, cp::FT)::FT
  (T_w - T_a) * ρ_a * cp * gH
end
