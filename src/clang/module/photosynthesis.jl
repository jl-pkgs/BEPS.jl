function photosynthesis(temp_leaf_p, rad_leaf, e_air, g_lb_w, vc_opt,
  f_soilwater, b_h2o, m_h2o, cii, temp_leaf_c, LH_leaf,
  Gs_w::TypeRef, aphoto::TypeRef, ci::TypeRef)

  ccall((:photosynthesis, libbeps), Cvoid,
    (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble,
      Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
    temp_leaf_p, rad_leaf, e_air, g_lb_w, vc_opt, f_soilwater, b_h2o, m_h2o, cii, temp_leaf_c, LH_leaf,
    Gs_w, aphoto, ci)
end


function photosynthesis(Tc_old::Leaf, R::Leaf, Ci_old::Leaf, leleaf::Leaf,
  temp_air, e_a10, f_soilwater, b_h2o, m_h2o,
  Gb_o::Cdouble, Gb_u::Cdouble, Vcmax_sunlit::Cdouble, Vcmax_shaded::Cdouble,
  Gs_new::LeafRef, Ac::LeafRef, Ci_new::LeafRef)

  photosynthesis(Tc_old.o_sunlit, R.o_sunlit, e_a10, Gb_o, Vcmax_sunlit, f_soilwater, b_h2o, m_h2o,
    Ci_old.o_sunlit,
    temp_air, leleaf.o_sunlit,
    Gs_new.o_sunlit, Ac.o_sunlit, Ci_new.o_sunlit)
  photosynthesis(Tc_old.o_shaded, R.o_shaded, e_a10, Gb_o, Vcmax_shaded, f_soilwater, b_h2o, m_h2o,
    Ci_old.o_shaded,
    temp_air, leleaf.o_shaded,
    Gs_new.o_shaded, Ac.o_shaded, Ci_new.o_shaded)
  photosynthesis(Tc_old.u_sunlit, R.u_sunlit, e_a10, Gb_u, Vcmax_sunlit, f_soilwater, b_h2o, m_h2o,
    Ci_old.u_sunlit,
    temp_air, leleaf.u_sunlit,
    Gs_new.u_sunlit, Ac.u_sunlit, Ci_new.u_sunlit)
  photosynthesis(Tc_old.u_shaded, R.u_shaded, e_a10, Gb_u, Vcmax_shaded, f_soilwater, b_h2o, m_h2o,
    Ci_old.u_shaded,
    temp_air, leleaf.u_shaded,
    Gs_new.u_shaded, Ac.u_shaded, Ci_new.u_shaded)
end
