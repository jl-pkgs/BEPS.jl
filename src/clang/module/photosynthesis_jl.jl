@fastmath function photosynthesis_jl(temp_leaf_p::Cdouble, Rsn_leaf::Cdouble, e_air::Cdouble, 
  g_lb_w::Cdouble, vc_opt::Cdouble,
  f_soilwater::Cdouble, b_h2o::Cdouble, m_h2o::Cdouble, 
  cii::Cdouble,
  T_leaf::Cdouble, LH_leaf::Cdouble)

  Gs_w::Cdouble = 0.0
  aphoto::Cdouble = 0.0
  ci::Cdouble = 0.0

  g_lb_c::Cdouble = 0.              # leaf laminar boundary layer condunctance to CO2 (mol m-2 s-1)
  rh_leaf::Cdouble = 0.             # relative humidity at leaf surface (0-1)
  temp_leaf_K::Cdouble = 0.         # leaf temperature (K)
  gs_co2_mole::Cdouble = 0.         # stomatal conductance to CO2 (mol m-2 s-1)
  gs_h2o_mole::Cdouble = 0.         # stomatal conductance to h2o (mol m-2 s-1)
  bc::Cdouble = 0.                  # temporary variable
  cs::Cdouble = 0.                  # CO2 concentration at leaf surface (ppm)
  b_co2::Cdouble = 0.      # the intercept term in BWB model (mol CO2 m-2 s-1): b_h2o/1.6
  m_co2::Cdouble = 0.      # the slope in BWB model: m_h2o/1.6
  gammac::Cdouble = 0.     # CO2 compensation point (ppm)
  jmopt::Cdouble = 0.      # the maximum potential electron transport rate at 25 deg C (umol m-2 s-1)
  jmax::Cdouble = 0.       # the maximum potential electron transport rate (umol m-2 s-1)
  vcmax::Cdouble = 0.      # the maximum velocities of carboxylation of Rubisco (umol m-2 s-1)
  km_co2::Cdouble = 0.     # Michaelis-Menten constant for CO2 (µmol mol-1)
  km_o2::Cdouble = 0.      # Michaelis-Menten constant for O2 (mmol mol-1)
  tau::Cdouble = 0.        # the specifity of Rubisco for CO2 compared with O2
  resp_ld::Cdouble = 0.    # leaf dark respiration (umol m-2 s-1)
  resp_ld25::Cdouble = 0.  # leaf dark respiration at 25 deg C (umol m-2 s-1)
  j_photon::Cdouble = 0.  # the flux of electrons through the thylakoid membrane (umol m-2 s-1)
  alpha_ps::Cdouble = 0.
  beta_ps::Cdouble = 0.
  gamma_ps::Cdouble = 0.
  theta_ps::Cdouble = 0.
  denom::Cdouble = 0.
  p_cubic::Cdouble = 0.
  q_cubic::Cdouble = 0.
  r_cubic::Cdouble = 0.
  Qroot::Cdouble = 0.
  Rroot::Cdouble = 0.

  # double root1, root2, root3;
  # double root_min = 0, root_max = 0, root_mid = 0;
  ang_L::Cdouble = 0.
  j_sucrose::Cdouble = 0.        # net photosynthesis rate limited by sucrose synthesis (umol m-2 s-1)
  # double wc, wj, psguess;  # gross photosynthesis rate limited by light (umol m-2 s-1)
  # double Aquad, Bquad, Cquad;
  # double b_ps, a_ps, e_ps, d_ps;
  product::Cdouble = 0.
  ps_1::Cdouble = 0.
  delta_1::Cdouble = 0.
  r3q::Cdouble = 0.
  tprime25::Cdouble = 0.

  ca::Cdouble = CO2_air                    # atmospheric co2 concentration (ppm)
  iphoton::Cdouble = 4.55 * 0.5 * Rsn_leaf # incident photosynthetic photon flux density (PPFD) umol m-2 s-1
  if (2 * iphoton < 1)
    iphoton = 0.
  end

  temp_leaf_K = T_leaf + 273.13

  fact_latent = LAMBDA(temp_leaf_p)
  bound_vapor = 1.0 / g_lb_w
  # air_pres::Cdouble = 101.325  # air pressure (kPa)
  #	g_lb_c = (g_lb_w/1.6)*air_pres/(temp_leaf_K*rugc); # (mol m-2 s-1)

  press_bars = 1.013
  pstat273 = 0.022624 / (273.16 * press_bars)

  T_Kelvin = T_leaf + 273.13
  rhova_g = e_air * 2165 / T_Kelvin  # absolute humidity, g m-3
  rhova_kg = rhova_g / 1000.0         # absolute humidity, kg m-3

  g_lb_c = 1.0 / (1.0 / g_lb_w * 1.6 * temp_leaf_K * (pstat273))

  m_co2 = m_h2o / 1.6
  b_co2 = b_h2o / 1.6

  rh_leaf = SFC_VPD(temp_leaf_K, LH_leaf, fact_latent, bound_vapor, rhova_kg)
  tprime25 = temp_leaf_K - tk_25  # temperature difference

  #/*****  Use Arrhenius Eq. to compute KC and km_o2  *****/
  km_co2 = TEMP_FUNC(kc25, ekc, tprime25, tk_25, temp_leaf_K)
  km_o2 = TEMP_FUNC(ko25, eko, tprime25, tk_25, temp_leaf_K)
  tau = TEMP_FUNC(tau25, ektau, tprime25, tk_25, temp_leaf_K)

  bc = km_co2 * (1.0 + o2 / km_o2)

  gammac = 0.5 * o2 / tau * 1000.0  # umol mol-1

  resp_ld25 = vc_opt * 0.004657

  # Bin Chen: check this later. reduce respiration by 40% in light according to Amthor
  if (2.0 * iphoton > 10)
    resp_ld25 *= 0.4
  end
  resp_ld = TEMP_FUNC(resp_ld25, erd, tprime25, tk_25, temp_leaf_K)

  #	jmopt = 29.1 + 1.64*vc_opt;
  jmopt = 2.39 * vc_opt - 14.2

  jmax = TBOLTZ(jmopt, ejm, toptjm, temp_leaf_K)    # Apply temperature correction to JMAX
  vcmax = TBOLTZ(vc_opt, evc, toptvc, temp_leaf_K)  # Apply temperature correction to vcmax

  # /*
  #  * APHOTO = PG - resp_ld, net photosynthesis is the difference
  #  * between gross photosynthesis and dark respiration. Note
  #  * photorespiration is already factored into PG.
  #  * **********************************************************
  #  *
  #  * Gs from Ball-Berry is for water vapor.  It must be divided
  #  * by the ratio of the molecular diffusivities to be valid for A
  #  */
  alpha_ps = 1.0 + (b_co2 / g_lb_c) - m_co2 * rh_leaf * f_soilwater
  beta_ps = ca * (g_lb_c * m_co2 * rh_leaf * f_soilwater - 2.0 * b_co2 - g_lb_c)
  gamma_ps = ca * ca * g_lb_c * b_co2
  theta_ps = g_lb_c * m_co2 * rh_leaf * f_soilwater - b_co2

  # /*
  #  * Test for the minimum of Wc and Wj.  Both have the form:
  #  *
  #  * W = (a ci - ad)/(e ci + b)
  #  *
  #  * after the minimum is chosen set a, b, e and d for the cubic solution.
  #  *
  #  * estimate of J according to Farquhar and von Cammerer (1981)
  #  */
  # /*if (jmax > 0)
  #     j_photon = qalpha * iphoton / sqrt(1. +(qalpha2 * iphoton * iphoton / (jmax * jmax)));
  # else
  #     j_photon = 0;*/
  # J photon from Harley
  j_photon = jmax * iphoton / (iphoton + 2.1 * jmax)

  # initial guess of intercellular CO2 concentration to estimate Wc and Wj:
  wj = j_photon * (cii - gammac) / (4.0 * cii + 8.0 * gammac)
  wc = vcmax * (cii - gammac) / (cii + bc)

  if (wj < wc)
    # for Harley and Farquhar type model for Wj
    psguess = wj
    a_ps = j_photon
    b_ps = 8.0 * gammac
    e_ps = 4.0
    d_ps = gammac
  else
    psguess = wc
    a_ps = vcmax
    b_ps = bc
    e_ps = 1.0
    d_ps = gammac
  end

  # if wj or wc are less than resp_ld then A would probably be less than zero.  This would yield a
  # negative stomatal conductance.  In this case, assume gs equals the cuticular value. This
  # assumptions yields a quadratic rather than cubic solution for A
  if (wj <= resp_ld || wc <= resp_ld)
    @goto quad
  end

  # cubic solution:
  #  A^3 + p A^2 + q A + r = 0
  denom = e_ps * alpha_ps

  p_cubic = (e_ps * beta_ps + b_ps * theta_ps - a_ps * alpha_ps + e_ps * resp_ld * alpha_ps)
  p_cubic /= denom

  q_cubic = (e_ps * gamma_ps + (b_ps * gamma_ps / ca) - a_ps * beta_ps + a_ps * d_ps * theta_ps +
             e_ps * resp_ld * beta_ps + resp_ld * b_ps * theta_ps)
  q_cubic /= denom

  r_cubic = -a_ps * gamma_ps + a_ps * d_ps * gamma_ps / ca +
            e_ps * resp_ld * gamma_ps + resp_ld * b_ps * gamma_ps / ca
  r_cubic /= denom

  # Use solution from Numerical Recipes from Press
  Qroot = (p_cubic * p_cubic - 3.0 * q_cubic) / 9.0
  Rroot = (2.0 * p_cubic * p_cubic * p_cubic - 9.0 * p_cubic * q_cubic + 27.0 * r_cubic) / 54.0

  (Qroot < 0) && (@goto quad)

  # @show Qroot
  r3q = Rroot / sqrt(Qroot * Qroot * Qroot)
  r3q = clamp(r3q, -1, 1) #  by G. Mo
  ang_L = acos(r3q)

  root1 = -2.0 * sqrt(Qroot) * cos(ang_L / 3.0) - p_cubic / 3.0  # real roots
  root2 = -2.0 * sqrt(Qroot) * cos((ang_L + PI2) / 3.0) - p_cubic / 3.0
  root3 = -2.0 * sqrt(Qroot) * cos((ang_L - PI2) / 3.0) - p_cubic / 3.0

  # Here A = x - p / 3, allowing the cubic expression to be expressed
  # as: x^3 + ax + b = 0
  # rank roots #1, #2 and #3 according to the minimum, intermediate and maximum value
  aphoto = findroot(root1, root2, root3)

  # also test for sucrose limitation of photosynthesis, as suggested by Collatz.  Js=Vmax/2
  j_sucrose = vcmax / 2.0 - resp_ld

  if (j_sucrose < aphoto)
    aphoto = j_sucrose
  end
  # Stomatal conductance for water vapor

  # Forests are hypostomatous.
  # Hence, we don't divide the total resistance
  # by 2 since transfer is going on only one side of a leaf.

  # if A < 0 then gs should go to cuticular value and recalculate A
  # using quadratic solution
  if aphoto <= 0.0
    @goto quad
  else
    @goto OUTDAT
  end
  # if aphoto < 0  set stomatal conductance to cuticle value
  # /*
  #  * a quadratic solution of A is derived if gs=b, but a cubic form occur
  #  * if gs = ax + b.  Use quadratic case when A <=0
  #  *
  #  * Bin Chen:
  #  * r_tot = 1.0/b_co2 + 1.0/g_lb_c; # total resistance to CO2 (m2 s mol-1)
  #  * denom = g_lb_c * b_co2;
  #  * Aquad = r_tot * e_ps;
  #  * Bquad = (e_ps*resp_ld + a_ps)*r_tot - b_ps - e_ps*ca;
  #  * Cquad = a_ps*(ca-d_ps) - resp_ld*(e_ps*ca+b_ps);
  #  */
  # original version
  @label quad
  ps_1 = ca * g_lb_c * b_co2
  delta_1 = b_co2 + g_lb_c
  denom = g_lb_c * b_co2

  Aquad = delta_1 * e_ps
  Bquad = -ps_1 * e_ps - a_ps * delta_1 + e_ps * resp_ld * delta_1 - b_ps * denom
  Cquad = a_ps * ps_1 - a_ps * d_ps * denom - e_ps * resp_ld * ps_1 - resp_ld * b_ps * denom

  product = Bquad * Bquad - 4.0 * Aquad * Cquad
  if (product >= 0)
    #	*aphoto = (-Bquad + sqrt(product)) / (2.0 * Aquad);
    aphoto = (-Bquad - sqrt(product)) / (2.0 * Aquad)
  end


  @label OUTDAT
  aphoto = max(0, aphoto)

  cs = ca - aphoto / g_lb_c

  gs_h2o_mole = (f_soilwater * m_h2o * rh_leaf * aphoto / cs) + b_h2o  # mol m-2 s-1
  gs_co2_mole = gs_h2o_mole / 1.6

  ci = cs - aphoto / gs_co2_mole
  Gs_w = gs_h2o_mole * temp_leaf_K * (pstat273)  # m s-1

  return Gs_w, aphoto, ci
end








@fastmath function SFC_VPD(temp_leaf_K::Float64, leleafpt::Float64,
  fact_latent::Float64, bound_vapor::Float64, rhova_kg::Float64)::Float64

  es_leaf = ES(temp_leaf_K)
  rhov_sfc = (leleafpt / (fact_latent)) * bound_vapor + rhova_kg # kg m-3
  e_sfc = rhov_sfc * temp_leaf_K / 0.2165 # mb
  vpd_sfc = es_leaf - e_sfc # mb
  rhum_leaf = 1.0 - vpd_sfc / es_leaf # 0 to 1.0
  return rhum_leaf
end

@fastmath function TEMP_FUNC(rate::Float64, eact::Float64, tprime::Float64, tref::Float64, t_lk::Float64)::Float64
  rate * exp(tprime * eact / (tref * rugc * t_lk))
end

# 
function LAMBDA(tak::Float64)::Float64
  y = 3149000.0 - 2370.0 * tak
  # add heat of fusion for melting ice
  if tak < 273.0
    y += 333.0
  end
  return y
end

# Function to calculate saturation vapor pressure function in mb
@fastmath function ES(t::Float64)::Float64
  if t > 0.0
    y1::Float64 = 54.8781919 - 6790.4985 / t - 5.02808 * log(t)
    y = exp(y1)
  else
    println("bad es calc")
    y = 0.0
  end
  return y
end

# Maxwell-Boltzmann temperature distribution for photosynthesis
@fastmath function TBOLTZ(rate::Float64, eakin::Float64, topt::Float64, tl::Float64)::Float64
  hkin::Float64 = 200000.0  # enthalpy term, J mol-1

  dtlopt::Float64 = tl - topt
  prodt::Float64 = rugc * topt * tl
  numm::Float64 = hkin * exp(eakin * dtlopt / prodt)
  denom::Float64 = hkin - eakin * (1.0 - exp(hkin * dtlopt / prodt))

  return rate * numm / denom
end


function findroot(root1::Float64, root2::Float64, root3::Float64)::Float64
  root_min::Float64 = 0.0
  root_mid::Float64 = 0.0
  root_max::Float64 = 0.0
  aphoto::Float64 = 0.0

  if root1 <= root2 && root1 <= root3
    root_min = root1
    if root2 <= root3
      root_mid = root2
      root_max = root3
    else
      root_mid = root3
      root_max = root2
    end
  end

  if root2 <= root1 && root2 <= root3
    root_min = root2
    if root1 <= root3
      root_mid = root1
      root_max = root3
    else
      root_mid = root3
      root_max = root1
    end
  end

  if root3 <= root1 && root3 <= root2
    root_min = root3
    if root1 < root2
      root_mid = root1
      root_max = root2
    else
      root_mid = root2
      root_max = root1
    end
  end

  # find out where roots plop down relative to the x-y axis
  if root_min > 0 && root_mid > 0 && root_max > 0
    aphoto = root_min
  end

  if root_min < 0 && root_mid < 0 && root_max > 0
    aphoto = root_max
  end

  if root_min < 0 && root_mid > 0 && root_max > 0
    aphoto = root_mid
  end
  # root_min, root_mid, root_max, 
  aphoto
end
