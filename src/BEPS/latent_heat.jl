function latent_heat!(leleaf::Leaf, Gw::Leaf, VPD, slope, Tc_old::Leaf, Tair, rho_a, Cp_ca, gamma)
  leleaf.o_sunlit = Gw.o_sunlit * (VPD + slope * (Tc_old.o_sunlit - Tair)) * rho_a * Cp_ca / gamma
  leleaf.o_shaded = Gw.o_shaded * (VPD + slope * (Tc_old.o_shaded - Tair)) * rho_a * Cp_ca / gamma
  leleaf.u_sunlit = Gw.u_sunlit * (VPD + slope * (Tc_old.u_sunlit - Tair)) * rho_a * Cp_ca / gamma
  leleaf.u_shaded = Gw.u_shaded * (VPD + slope * (Tc_old.u_shaded - Tair)) * rho_a * Cp_ca / gamma
end
