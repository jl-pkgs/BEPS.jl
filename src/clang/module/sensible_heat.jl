function sensible_heat(tempL::Leaf, 
  temp_g::Cdouble, temp_air::Cdouble, rh_air::Cdouble,
  Gheat::Leaf, Gheat_g::Cdouble, LAI::Leaf)

  SH_o = init_dbl()
  SH_u = init_dbl()
  SH_g = init_dbl()

  ccall((:sensible_heat, libbeps), Cvoid,
    (Leaf, Cdouble, Cdouble, Cdouble, Leaf, Cdouble, Leaf,
      Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
    tempL, temp_g, temp_air, rh_air, Gheat, Gheat_g, LAI,
    SH_o, SH_u, SH_g)
  
  SH_o[], SH_u[], SH_g[]
end


function sensible_heat_jl(T_leaf::Leaf, T_ground::FT, T_air::FT, rh_air::FT,
  Gheat::Leaf, Gheat_g::FT, lai::Leaf)

  SH = Leaf()
  # meteo_pack_output = meteo_pack(T_air, rh_air)
  q = RH2q(T_air, rh_air)
  cp_air = cal_cp(q)
  rho_air0 = rho_a

  SH.o_sunlit = (T_leaf.o_sunlit - T_air) * rho_air0 * cp_air * Gheat.o_sunlit
  SH.o_shaded = (T_leaf.o_shaded - T_air) * rho_air0 * cp_air * Gheat.o_shaded
  SH.u_sunlit = (T_leaf.u_sunlit - T_air) * rho_air0 * cp_air * Gheat.u_sunlit
  SH.u_shaded = (T_leaf.u_shaded - T_air) * rho_air0 * cp_air * Gheat.u_shaded

  SH_o = SH.o_sunlit * lai.o_sunlit + SH.o_shaded * lai.o_shaded
  SH_u = SH.u_sunlit * lai.u_sunlit + SH.u_shaded * lai.u_shaded

  SH_o = max(-200, SH_o)
  SH_u = max(-200, SH_u)
  SH_g = (T_ground - T_air) * rho_air0 * cp_air * Gheat_g

  SH_o, SH_u, SH_g
end

# function cal_SH(T_leaf, T_air, Gheat, rho_air0, cp_air)
#   (T_leaf - T_air) * rho_air0 * cp_air * Gheat
# end
