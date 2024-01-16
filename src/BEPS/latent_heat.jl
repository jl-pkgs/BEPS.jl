function latent_heat!(leleaf::Leaf, Gw::Leaf, VPD_air, slope, Tc_old::Leaf, temp_air, rho_a, Cp_ca, psychrometer)
  leleaf.o_sunlit = Gw.o_sunlit * (VPD_air + slope * (Tc_old.o_sunlit - temp_air)) * rho_a * Cp_ca / psychrometer
  leleaf.o_shaded = Gw.o_shaded * (VPD_air + slope * (Tc_old.o_shaded - temp_air)) * rho_a * Cp_ca / psychrometer
  leleaf.u_sunlit = Gw.u_sunlit * (VPD_air + slope * (Tc_old.u_sunlit - temp_air)) * rho_a * Cp_ca / psychrometer
  leleaf.u_shaded = Gw.u_shaded * (VPD_air + slope * (Tc_old.u_shaded - temp_air)) * rho_a * Cp_ca / psychrometer
end
