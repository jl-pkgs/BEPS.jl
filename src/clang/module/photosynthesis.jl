include("photosynthesis_jl.jl")


function photosynthesis_c(temp_leaf_p::Cdouble, rad_leaf::Cdouble, e_air::Cdouble, 
  g_lb_w::Cdouble, vc_opt::Cdouble,
  f_soilwater::Cdouble, b_h2o::Cdouble, m_h2o::Cdouble, 
  cii::Cdouble, temp_leaf_c::Cdouble, LH_leaf::Cdouble)

  Gs_w = init_dbl()
  aphoto = init_dbl()
  ci = init_dbl()

  ccall((:photosynthesis, libbeps), Cvoid,
    (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble,
      Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
    temp_leaf_p, rad_leaf, e_air, g_lb_w, vc_opt, f_soilwater, b_h2o, m_h2o, cii, temp_leaf_c, LH_leaf,
    Gs_w, aphoto, ci)

  Gs_w[], aphoto[], ci[]
end


function photosynthesis(Tc_old::Leaf, R::Leaf, Ci_old::Leaf, leleaf::Leaf,
  Ta::Cdouble, ea::Cdouble, f_soilwater::Cdouble, b_h2o::Cdouble, m_h2o::Cdouble,
  Gb_o::Cdouble, Gb_u::Cdouble, Vcmax_sunlit::Cdouble, Vcmax_shaded::Cdouble,
  # output
  Gs_new::Leaf, Ac::Leaf, Ci_new::Leaf; version = "c")
  
  if version == "c"
    fun = photosynthesis_c
  elseif version == "julia"
    fun = photosynthesis_jl
  end
  
  Gs_new.o_sunlit, Ac.o_sunlit, Ci_new.o_sunlit =
    fun(Tc_old.o_sunlit, R.o_sunlit, ea, Gb_o, Vcmax_sunlit, f_soilwater, b_h2o, m_h2o,
      Ci_old.o_sunlit,
      Ta, leleaf.o_sunlit)

  Gs_new.o_shaded, Ac.o_shaded, Ci_new.o_shaded =
    fun(Tc_old.o_shaded, R.o_shaded, ea, Gb_o, Vcmax_shaded, f_soilwater, b_h2o, m_h2o,
      Ci_old.o_shaded,
      Ta, leleaf.o_shaded)

  Gs_new.u_sunlit, Ac.u_sunlit, Ci_new.u_sunlit =
    fun(Tc_old.u_sunlit, R.u_sunlit, ea, Gb_u, Vcmax_sunlit, f_soilwater, b_h2o, m_h2o,
      Ci_old.u_sunlit,
      Ta, leleaf.u_sunlit)

  Gs_new.u_shaded, Ac.u_shaded, Ci_new.u_shaded =
    fun(Tc_old.u_shaded, R.u_shaded, ea, Gb_u, Vcmax_shaded, f_soilwater, b_h2o, m_h2o,
      Ci_old.u_shaded,
      Ta, leleaf.u_shaded)
end
