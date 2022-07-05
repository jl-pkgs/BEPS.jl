#/ @brief the inter-module function between main program and modules
#/ @param  jday       day of year
#/ @param  rstep      hour of day
#/ @param  lai        leaf area index
#/ @param  clumping   clumping index
#/ @param  parameter  parameter array according to land cover types
#/ @param  meteo      meteorological data
#/ @param  CosZs      cosine of solar zenith angle
#/ @param  var_o      temporary variables array of last time step
#/ @param  var_n      temporary variables array of this time step
#/ @param  soilp      soil coefficients according to land cover types and soil textures
#/ @param  mid_res    results struct
#/ @return void

# kPa deg-1
cal_slope(Ta::Float64) = 2503.0 / pow((Ta + 237.3), 2) * exp(17.27 * Ta / (Ta + 237.3))

# kPa
cal_es(Ta::Float64) = 0.61078 * exp(17.3 * Ta / (237.3 + Ta))

function inter_prg(jday::Int, rstep::Int, lai::Float64, clumping::Float64, parameter::Vector{T}, meteo::ClimateData,
  CosZs::Float64,
  var_o::Vector{T}, var_n::Vector{T}, soilp::Soil, mid_res::Results) where {T}
  # /*****  define parameters and arrays  *****/
  # int num,kkk;
  #   int landcover;
  #   double lai_o,lai_u;
  # double stem_o,stem_u;

  #   double d_soil[layer+1];
  #   # double Zsp;                               # the depth of snow on the surface
  #   double Zp,Zp1=0,Zp2=0;                    # depth of pounded water on the surface
  #   double height_wind_sp;                    # height of the Va measured for calculation of L

  #   double Qhc_o[120] = {0}, Qhc_u[120] = {0}, Qhg[120] = {0};  # The sensible heat flux from canopy and ground
  #   double G[layer + 2][120] = {0};                             # the heat flux into the canopy of over story --in W/m^2

  #   double Wcl_o[120] = {0}, Wcs_o[120] = {0};  # the masses od rain and snow on the canopy
  #   double Xcl_o[120] = {0}, Xcs_o[120] = {0};  # the fraction of canopy covered by liquid water and snow
  #   double Wcl_u[120] = {0}, Wcs_u[120] = {0};  # the masses of rain and snow on the canopy
  #   double Xcl_u[120] = {0}, Xcs_u[120] = {0};  # the fraction of canopy covered by liquid water and snow

  #   double r_rain_g[120] = {0};               # the rainfall rate, on ground surface  m/s
  #   double rho_snow[120] = {0};               # density of snow
  #   double alpha_v_sw[120] = {0}, alpha_n_sw[120] = {0};  # albedo of snow
  #   double Wg_snow[120] = {0};                            # the amount of snow on the ground
  #   double Xg_snow[120] = {0};                            # the fraction of the ground surface covered by snow
  #   double Ac_snow_u[120] = {0};                          # the areas of canopy covered in snow, o--overstory and  u-- understory
  #   double Ac_snow_o[120] = {0};

  #   double Ts0[120] = {0}, Tsn0[120] = {0}, Tsm0[120] = {0}, Tsn1[120] = {0}, Tsn2[120] = {0};  # surface temperature
  #   double Tc_u[120] = {0};                                             # the effective canopy temperature in K
  #   double Tm[layer + 2][120] = {0};                                    # Tb[layer+2][120],soil temperature at the bottom and the middle of each layer*/

  #   double lambda[layer + 2] = {0};  # thermal conductivity of each soil layer;*/
  #   double Cs[layer + 2][120] = {0};  # the soil volumetric heat capacity of each soil layer, j+kkk/m^3/K*/

  #   double Ta;                 # air temperature */
  #   double precip, rh_air, wind_sp;  # precipitation in meters, relative humidity (%), wind speed in m/s */
  #   double temp_grd;                 # ground temperature */

  #   double Eil_o[120] = {0}, EiS_o[120] = {0};  # the evaporation rate of intercepted water of overstory--in kg/m^2/s; l-- water; S-snow
  #   double Eil_u[120] = {0}, EiS_u[120] = {0};  # the evaporation rate of intercepted water of overstory--in kg/m^2/s; l-- water; S-snow of intercepted water--in kg/m^2/s
  #   double Trans_o[120] = {0}, Trans_u[120] = {0};  # transpiration from overstory and understory canopies
  #   double Evap_soil[120] = {0};                    # evaporation from soil
  #   double Evap_SW[120] = {0};                      # evaporation from water pond
  #   double Evap_SS[120] = {0};                      # evaporation from snow pack

  #   double lambda_snow[120] = {0};  # the effective thermal conductivity of snow --in m^2/s
  #   double e_a10;             # the vapour partial pressure of water --in kPa(1mb=100Pa=0.1kpa)

  #   double Lv_liquid;                  # the latent heat of evaporation from liquid at air temperature=Ta, in j+kkk/kg
  #   double Lv_solid = 2.83 * 1000000;  # the latent heat of evaporation from solid (snow/ice) at air temperature=Ta, in j+kkk/kg

  #   double Ks;  # KsCal,KsMea[120] instantaneous total short wave radiation (Global radiation)
  #   double alpha_sat,alpha_dry;
  #   double alpha_v_o,alpha_v_u;               # visible albedo of overstory,  o--overstory, u--understory;
  #   double alpha_n_o,alpha_n_u;               # near-infrared albedo of overstory,o--overstory, u--understory;
  #   double alpha_g;                           # the all-wave ground surface albedo
  #   double alpha_v_g, alpha_n_g;              # the ground surface albedo for visible range and for near-infrared respectively
  #   double Cp_ca;                             # specific heat of moist air above the canopy  in j+kkk/k/G
  #   double ra_o,ra_u,ra_g;                    # the aerodynamic resistance of overstory, understory and ground surface   in s/m

  # double q_ca;                              # the actual canopy stomatal resistance  --in s/m
  # double radiation_o, radiation_u, radiation_g;

  # double ip=0;                              # the cumulative infiltration at the time of ponding   --in m/s
  # double Inf=0;
  # double zr=0.8;
  # double Cpd=1004.65;
  # double rho_w=1025.0;	                  # density of water

  # double Cs_o_sunlit_old,Cs_o_shaded_old,Cs_u_sunlit_old,Cs_u_shaded_old; # CO2 concentration on the surfaces of leaves
  # double Tc_o_sunlit_old,Tc_o_shaded_old,Tc_u_sunlit_old,Tc_u_shaded_old; # the effective canopy temperature in K

  # double Gs_o_sunlit_new,Gs_o_shaded_new,Gs_u_sunlit_new,Gs_u_shaded_new; # stomatal conductance of the big leaf for water
  # double Gs_o_sunlit_old,Gs_o_shaded_old,Gs_u_sunlit_old,Gs_u_shaded_old; # stomatal conductance of the big leaf for water

  # double Ac_o_sunlit,Ac_o_shaded,Ac_u_sunlit,Ac_u_shaded;                 # net photosynthesis rate
  # double Cs_o_sunlit_new,Cs_o_shaded_new,Cs_u_sunlit_new,Cs_u_shaded_new; # CO2 concentration on the surfaces of leaves
  # double Ci_o_sunlit_new,Ci_o_shaded_new,Ci_u_sunlit_new,Ci_u_shaded_new; # intercellular CO2 concentration pn the leaf
  # double Ci_o_sunlit_old,Ci_o_shaded_old,Ci_u_sunlit_old,Ci_u_shaded_old; # intercellular CO2 concentration pn the leaf
  # double Cc_o_sunlit_new,Cc_o_shaded_new,Cc_u_sunlit_new,Cc_u_shaded_new; # CO2 concentration in the chloroplast

  # double Tc_o_sunlit_new,Tc_o_shaded_new,Tc_u_sunlit_new,Tc_u_shaded_new; # the effective canopy temperature in K
  # double f_soilwater;                                                     # an empirical parameter describing the relative availability of soil water for plants
  # double Gw_o_sunlit,Gw_o_shaded,Gw_u_sunlit,Gw_u_shaded;                 # the total conductance for water from the intercellular space of the leaves to the reference height above the canopy
  # double Gc_o_sunlit,Gc_o_shaded,Gc_u_sunlit,Gc_u_shaded;                 # the total conductance for CO2 from the intercellular space of the leaves to the reference height above the canopy
  # double Gww_o_sunlit,Gww_o_shaded,Gww_u_sunlit,Gww_u_shaded;             # the total conductance for water from the surface of the leaves to the reference height above the canopy
  # double Gh_o_sunlit,Gh_o_shaded,Gh_u_sunlit,Gh_u_shaded;                 # total conductance for heat transfer from the leaf surface to the reference height above the canopy

  # double psychrometer;
  # double R_o_sunlit,R_o_shaded,R_u_sunlit,R_u_shaded;                     # solar radiation absorbed by sunlit, shaded leaves

  # double Tco, Tcu,slope;
  # double H_o_sunlit,H_o_shaded;                                           # sensible heat flux from leaves
  # double LAI_o_sunlit,LAI_o_shaded,LAI_u_sunlit,LAI_u_shaded;             # stem_o,stem_u;
  # double LAIo_sunlit,LAIo_shaded,LAIu_sunlit,LAIu_shaded;
  # double radiation_o_sun, radiation_o_shaded;                             # net radiation of leaves
  # double radiation_u_sun, radiation_u_shaded;
  # double GPP_o_sunlit,GPP_o_shaded,GPP_u_sunlit,GPP_u_shaded;

  # double VPS_air;
  # double GH_o,G_o_a, G_o_b,G_u_a, G_u_b;
  # double canopyh_o,canopyh_u;

  # double VPD_air;             # Vapor pressure deficit term

  # double mass_water_g;
  # double percentArea_snow_o, percentArea_snow_u;
  # double Gheat_g;

  # double b_h2o;               # the intercept term in BWB model (mol H2O m-2 s-1)
  # double m_h2o;               # the slope in BWB model

  #   # leaf latent heat flux (mol/m2/s)
  # double leleaf_o_sunlit;
  # double leleaf_o_shaded;
  # double leleaf_u_sunlit;
  # double leleaf_u_shaded;

  # # parameters for Vcmax-Nitrogen calculations
  Kn = 0.3  #0.713/2.4
  G_theta = 0.5
  # double K,Vcmax0,Vcmax_sunlit,Vcmax_shaded,expr1,expr2,expr3;
  # double slope_Vcmax_N, leaf_N,Jmax_sunlit,Jmax_shaded;
  psychrometer = 0.066

  alpha_sat = parameter[24+1]       # albedo of saturated/dry soil for module rainfall 1
  alpha_dry = parameter[25+1]       # the albedo of dry soil
  canopyh_o = parameter[29+1]     # to be used for module aerodynamic_conductance
  canopyh_u = parameter[30+1]
  height_wind_sp = parameter[31+1] # the_height_to_measure_wind_speed, for module aerodynamic_conductance
  m_h2o = parameter[33+1]           # to be used for module photosynthesis
  b_h2o = parameter[34+1]

  # /*****  Vcmax-Nitrogen calculations，by G.Mo，Apr. 2011  *****/
  if (CosZs > 0) # day time
    K = G_theta * clumping / CosZs  # G_theta = 0.5 assuming a spherical leaf angle distribution
    Vcmax0 = parameter[36+1]
    expr1 = 1 - exp(-K * lai)
    expr2 = 1 - exp(-lai * (Kn + K))
    expr3 = 1 - exp(-Kn * lai)

    # Formulas based on Chen et al., 2012, GBC
    if (expr1 > 0)
      Vcmax_sunlit = Vcmax0 * parameter[47+1] * parameter[46+1] * K * expr2 / (Kn + K) / expr1
    else
      Vcmax_sunlit = Vcmax0
    end

    if (K > 0 && lai > expr1 / K)
      Vcmax_shaded = Vcmax0 * parameter[47+1] * parameter[46+1] * (expr3 / Kn - expr2 / (Kn + K)) / (lai - expr1 / K)
    else
      Vcmax_shaded = Vcmax0
    end
  end

  # /*****  LAI calculation module, by B. Chen  *****/
  lai_o = lai
  if (lai < 0.1)
    lai_o = 0.1
  end
  landcover = Int(parameter[4+1])

  if (landcover == 25 || landcover == 40)
    lai_u = 0.01
  else
    lai_u = 1.18 * exp(-0.99 * lai_o)
  end

  if (lai_u > lai_o)
    lai_u = 0.01
  end

  stem_o = parameter[8+1] * 0.2    # parameter[8].LAI max overstory
  stem_u = parameter[9+1] * 0.2    # parameter[9].LAI max understory

  # lai_calc module
  # separate lai into sunlit and shaded portions
  lai2(stem_o, stem_u, landcover, CosZs, lai_o, clumping, lai_u,
    LAIo_sunlit, LAIo_shaded, LAIu_sunlit, LAIu_shaded, # Ref
    LAI_o_sunlit, LAI_o_shaded, LAI_u_sunlit, LAI_u_shaded)

  # /*****  Initialization of this time step  *****/
  Ks = meteo.Srad
  rh_air = meteo.rh
  wind_sp = meteo.wind
  precip = meteo.rain / step  # precipitation in meters
  Ta = meteo.temp

  if (Ks <= 0)
    alpha_v_o = 0
    alpha_n_o = 0
    alpha_v_u = 0
    alpha_n_u = 0
  else
    alpha_v_o = parameter[22+1]
    alpha_n_o = parameter[23+1]
    alpha_v_u = parameter[22+1]
    alpha_n_u = parameter[23+1]
  end

  # Ground surface temperature
  Ts0[1] = clamp(var_o[3+1], Ta - 2.0, Ta + 2.0)
  Tsn0[1] = clamp(var_o[4+1], Ta - 2.0, Ta + 2.0)
  Tsm0[1] = clamp(var_o[5+1], Ta - 2.0, Ta + 2.0)
  Tsn1[1] = clamp(var_o[6+1], Ta - 2.0, Ta + 2.0)
  Tsn2[1] = clamp(var_o[7+1], Ta - 2.0, Ta + 2.0)

  # if ((Ts0[0]-Ta)>2.0) Ts0[0]=Ta+2.0;
  # if ((Ts0[0]-Ta)<-2.0) Ts0[0]=Ta-2.0;
  # Tsn0[0] =var_o[4] ;
  # if ((Tsn0[0]-Ta)>2.0) Tsn0[0]=Ta+2.0;
  # if ((Tsn0[0]-Ta)<-2.0) Tsn0[0]=Ta-2.0;
  # Tsm0[0] =var_o[5];
  # if ((Tsm0[0]-Ta)>2.0) Tsm0[0]=Ta+2.0;
  # if ((Tsm0[0]-Ta)<-2.0) Tsm0[0]=Ta-2.0;
  # Tsn1[0] =var_o[6] ;
  # if ((Tsn1[0]-Ta)>2.0) Tsn1[0]=Ta+2.0;
  # if ((Tsn1[0]-Ta)<-2.0) Tsn1[0]=Ta-2.0;
  # Tsn2[0] =var_o[7] ;
  # if ((Tsn2[0]-Ta)>2.0) Tsn2[0]=Ta+2.0;
  # if ((Tsn2[0]-Ta)<-2.0) Tsn2[0]=Ta-2.0;
  Qhc_o[1] = var_o[11+1]
  Wcl_o[1] = var_o[15+1]
  Wcs_o[1] = var_o[16+1]   # /* the mass of intercepted liquid water and snow, overstory */

  # the evaporation rate of rain and snow--in kg/m^2/s, understory
  Wcl_u[1] = var_o[18+1]
  Wcs_u[1] = var_o[19+1]   # /* the mass of intercepted liquid water and snow, overstory */
  Wg_snow[1] = var_o[20+1]  # /* thr fraction of ground surface covered in snow and snow mass */

  Zsp = soilp.Zsp
  Zp = soilp.Zp

  if (Zp < 0.001)
    Zp = 0
  end

  # /*****  Vcmax Jmax module by L. He  *****/
  #slope_Vcmax_N = parameter[47];
  #leaf_N = parameter[46];

  #Vcmax_Jmax(lai_o, clumping, Vcmax0,slope_Vcmax_N, leaf_N, CosZs, &Vcmax_sunlit, &Vcmax_shaded, &Jmax_sunlit, &Jmax_shaded);

  # temperatures of overstorey and understorey canopies
  Tc_o_sunlit_old = Ta - 0.5
  Tc_o_shaded_old = Ta - 0.5
  Tc_u_sunlit_old = Ta - 0.5
  Tc_u_shaded_old = Ta - 0.5

  # /*****  Ten time intervals in a hourly time step.6min or 360s per loop  ******/
  # for(kkk = 1;kkk <= kloop;kkk++)
  for kkk = 1:kloop
    # /*****  Snow pack stage 1 by X. Luo  *****/
    snowpack_stage1(Ta, precip, Wcs_o[kkk-1], Wcs_u[kkk-1], Wg_snow[kkk-1],
      &Wcs_o[kkk], &Wcs_u[kkk], &Wg_snow[kkk],
      lai_o, lai_u, clumping,
      &Ac_snow_o[kkk], &Ac_snow_u[kkk], &Xcs_o[kkk], &Xcs_u[kkk], &Xg_snow[kkk],
      &rho_snow[kkk], &Zsp, &alpha_v_sw[kkk], &alpha_n_sw[kkk])

    # /*****  Rain fall stage 1 by X. Luo  *****/
    rainfall_stage1(Ta, precip, Wcl_o[kkk-1], Wcl_u[kkk-1],
      lai_o, lai_u, clumping, &Wcl_o[kkk], &Wcl_u[kkk], &Xcl_o[kkk], &Xcl_u[kkk], &r_rain_g[kkk])

    # Old version
    # if(thetam[0][kkk-1]<soilp.theta_vwp[1]*0.5) alpha_g = alpha_dry;
    # else alpha_g = (thetam[0][kkk-1]-soilp.theta_vwp[1]*0.5)/(soilp.fei[1]-soilp.theta_vwp[1]*0.5) * (alpha_sat - alpha_dry) + alpha_dry;
    if (soilp.thetam_prev[1] < soilp.theta_vwp[1] * 0.5)
      alpha_g = alpha_dry
    else
      alpha_g = (soilp.thetam_prev[1] - soilp.theta_vwp[1] * 0.5) / (soilp.fei[1] - soilp.theta_vwp[1] * 0.5) * (alpha_sat - alpha_dry) + alpha_dry
    end

    alpha_v_g = 2.0 / 3.0 * alpha_g
    alpha_n_g = 4.0 / 3.0 * alpha_g


    # /*****  Soil water factor module by L. He  *****/
    soil_water_factor_v2(soilp)
    f_soilwater = soilp.f_soilwater

    if (f_soilwater > 1.0)
      f_soilwater = 1.0
    end  # to be used for module photosynthesis

    GH_o = Qhc_o[kkk-1]# to be used as the init. for module aerodynamic_conductance

    VPS_air = cal_es(Ta)  # to estimate saturated water vapor pressure in kpa
    e_a10 = VPS_air * rh_air / 100                      # to be used for module photosynthesis
    VPD_air = VPS_air - e_a10                              # water vapor deficit at the reference height

    q_ca = 0.622 * e_a10 / (101.35 - 0.378 * e_a10)      # in g/g, unitless
    Cp_ca = Cpd * (1 + 0.84 * q_ca)

    slope = cal_slope(Ta)

    Gs_o_sunlit_old = 1 / 200.0
    Ci_o_sunlit_old = 0.7 * CO2_air
    Gs_o_shaded_old = 1 / 200.0
    Ci_o_shaded_old = 0.7 * CO2_air
    Gs_u_sunlit_old = 1 / 300.0
    Ci_u_sunlit_old = 0.7 * CO2_air
    Gs_u_shaded_old = 1 / 300.0
    Ci_u_shaded_old = 0.7 * CO2_air

    percentArea_snow_o = Ac_snow_o[kkk] / lai_o / 2
    percentArea_snow_u = Ac_snow_u[kkk] / lai_u / 2

    temp_grd = Ta   # ground temperature substituted by air temperature

    num = 0
    while true # iteration for BWB equation until results converge
      num = num + 1

      # /***** Aerodynamic_conductance module by G.Mo  *****/

      aerodynamic_conductance(canopyh_o, canopyh_u, height_wind_sp, clumping, Ta, wind_sp, GH_o,
        lai_o + stem_o, lai_u + stem_u, &ra_o, &ra_u, &ra_g, &G_o_a, &G_o_b, &G_u_a, &G_u_b)


      Gh_o_sunlit = 1.0 / (1.0 / G_o_a + 0.5 / G_o_b)      # heat conductance of sunlit leaves of overstorey
      Gh_o_shaded = 1.0 / (1.0 / G_o_a + 0.5 / G_o_b)      # heat conductance of shaded leaves of overstorey
      Gh_u_sunlit = 1.0 / (1.0 / G_u_a + 0.5 / G_u_b)      # heat conductance of sunlit leaves of understorey
      Gh_u_shaded = 1.0 / (1.0 / G_u_a + 0.5 / G_u_b)      # heat conductance of shaded leaves of understorey

      Gww_o_sunlit = 1.0 / (1.0 / G_o_a + 1.0 / G_o_b + 100) # conductance for intercepted water of sunlit leaves of overstorey
      Gww_o_shaded = 1.0 / (1.0 / G_o_a + 1.0 / G_o_b + 100) # conductance for intercepted water of shaded leaves of overstorey
      Gww_u_sunlit = 1.0 / (1.0 / G_u_a + 1.0 / G_u_b + 100) # conductance for intercepted water of sunlit leaves of understorey
      Gww_u_shaded = 1.0 / (1.0 / G_u_a + 1.0 / G_u_b + 100) # conductance for intercepted water of shaded leaves of understorey

      # temperatures of overstorey and understorey canopies
      Tco = (Tc_o_sunlit_old * LAI_o_sunlit + Tc_o_shaded_old * LAI_o_shaded) / (LAI_o_sunlit + LAI_o_shaded)
      Tcu = (Tc_u_sunlit_old * LAI_u_sunlit + Tc_u_shaded_old * LAI_u_shaded) / (LAI_u_sunlit + LAI_u_shaded)


      # /*****  Net radiation at canopy and leaf level module by X.Luo  *****/

      netRadiation(Ks, CosZs, Tco, Tcu, temp_grd, lai_o, lai_u, lai_o + stem_o, lai_u + stem_u, LAI_o_sunlit, LAI_o_shaded, LAI_u_sunlit, LAI_u_shaded,
        clumping, Ta, rh_air, alpha_v_sw[kkk], alpha_n_sw[kkk], percentArea_snow_o, percentArea_snow_u,
        Xg_snow[kkk], alpha_v_o, alpha_n_o, alpha_v_u, alpha_n_u, alpha_v_g, alpha_n_g,
        &radiation_o, &radiation_u, &radiation_g, &radiation_o_sun, &radiation_o_shaded, &radiation_u_sun, &radiation_u_shaded,
        &R_o_sunlit, &R_o_shaded, &R_u_sunlit, &R_u_shaded)


      # /*****  Photosynthesis module by B. Chen  *****/

      # Four components: overstory sunlit and shaded, understory sunlit and shaded
      Gw_o_sunlit = 1.0 / (1.0 / G_o_a + 1.0 / G_o_b + 1.0 / Gs_o_sunlit_old) # conductance of sunlit leaves of overstorey for water
      Gw_o_shaded = 1.0 / (1.0 / G_o_a + 1.0 / G_o_b + 1.0 / Gs_o_shaded_old) # conductance of shaded leaves of overstorey for water
      Gw_u_sunlit = 1.0 / (1.0 / G_u_a + 1.0 / G_u_b + 1.0 / Gs_u_sunlit_old) # conductance of sunlit leaves of understorey for water
      Gw_u_shaded = 1.0 / (1.0 / G_u_a + 1.0 / G_u_b + 1.0 / Gs_u_shaded_old) # conductance of shaded leaves of understorey for water


      leleaf_o_sunlit = Gw_o_sunlit * (VPD_air + slope * (Tc_o_sunlit_old - Ta)) * rho_a * Cp_ca / psychrometer
      leleaf_o_shaded = Gw_o_shaded * (VPD_air + slope * (Tc_o_shaded_old - Ta)) * rho_a * Cp_ca / psychrometer
      leleaf_u_sunlit = Gw_u_sunlit * (VPD_air + slope * (Tc_u_sunlit_old - Ta)) * rho_a * Cp_ca / psychrometer
      leleaf_u_shaded = Gw_u_shaded * (VPD_air + slope * (Tc_u_shaded_old - Ta)) * rho_a * Cp_ca / psychrometer

      if (CosZs > 0)

        photosynthesis(Tc_o_sunlit_old, R_o_sunlit, e_a10, G_o_b, Vcmax_sunlit, f_soilwater, b_h2o, m_h2o, Ci_o_sunlit_old,
          Ta, leleaf_o_sunlit, &Gs_o_sunlit_new, &Ac_o_sunlit, &Ci_o_sunlit_new)
        photosynthesis(Tc_o_shaded_old, R_o_shaded, e_a10, G_o_b, Vcmax_shaded, f_soilwater, b_h2o, m_h2o, Ci_o_shaded_old,
          Ta, leleaf_o_shaded, &Gs_o_shaded_new, &Ac_o_shaded, &Ci_o_shaded_new)
        photosynthesis(Tc_u_sunlit_old, R_u_sunlit, e_a10, G_u_b, Vcmax_sunlit, f_soilwater, b_h2o, m_h2o, Ci_u_sunlit_old,
          Ta, leleaf_u_sunlit, &Gs_u_sunlit_new, &Ac_u_sunlit, &Ci_u_sunlit_new)
        photosynthesis(Tc_u_shaded_old, R_u_shaded, e_a10, G_u_b, Vcmax_shaded, f_soilwater, b_h2o, m_h2o, Ci_u_shaded_old,
          Ta, leleaf_u_shaded, &Gs_u_shaded_new, &Ac_u_shaded, &Ci_u_shaded_new)
      else
        Gs_o_sunlit_new = 0.0001
        Ac_o_sunlit = 0.0
        Ci_o_sunlit_new = CO2_air * 0.7
        Cs_o_sunlit_new = CO2_air
        Cc_o_sunlit_new = CO2_air * 0.7 * 0.8

        Gs_o_shaded_new = 0.0001
        Ac_o_shaded = 0.0
        Ci_o_shaded_new = CO2_air * 0.7
        Cs_o_shaded_new = CO2_air
        Cc_o_shaded_new = CO2_air * 0.7 * 0.8

        Gs_u_sunlit_new = 0.0001
        Ac_u_sunlit = 0.0
        Ci_u_sunlit_new = CO2_air * 0.7
        Cs_u_sunlit_new = CO2_air
        Cc_u_sunlit_new = CO2_air * 0.7 * 0.8

        Gs_u_shaded_new = 0.0001
        Ac_u_shaded = 0.0
        Ci_u_shaded_new = CO2_air * 0.7
        Cs_u_shaded_new = CO2_air
        Cc_u_shaded_new = CO2_air * 0.7 * 0.8
      end

      Ci_o_sunlit_old = Ci_o_sunlit_new
      Cs_o_sunlit_old = Cs_o_sunlit_new
      Gs_o_sunlit_old = Gs_o_sunlit_new # m/s
      Gw_o_sunlit = 1.0 / (1.0 / G_o_a + 1.0 / G_o_b + 1.0 / Gs_o_sunlit_new) #conductance of sunlit leaves of overstorey for water
      Gc_o_sunlit = 1.0 / (1.0 / G_o_a + 1.4 / G_o_b + 1.6 / Gs_o_sunlit_new) #conductance of sunlit leaves of overstorey for CO2

      Ci_o_shaded_old = Ci_o_shaded_new
      Cs_o_shaded_old = Cs_o_shaded_new
      Gs_o_shaded_old = Gs_o_shaded_new   # m/s
      Gw_o_shaded = 1.0 / (1.0 / G_o_a + 1.0 / G_o_b + 1.0 / Gs_o_shaded_new) #conductance of sunlit leaves of overstorey for water
      Gc_o_shaded = 1.0 / (1.0 / G_o_a + 1.4 / G_o_b + 1.6 / Gs_o_shaded_new) #conductance of sunlit leaves of overstorey for CO2


      Ci_u_sunlit_old = Ci_o_sunlit_new
      Cs_u_sunlit_old = Cs_u_sunlit_new
      Gs_u_sunlit_old = Gs_o_sunlit_new   # m/s
      Gw_u_sunlit = 1.0 / (1.0 / G_u_a + 1.0 / G_u_b + 1.0 / Gs_u_sunlit_new) #conductance of sunlit leaves of overstorey for water
      Gc_u_sunlit = 1.0 / (1.0 / G_u_a + 1.4 / G_u_b + 1.6 / Gs_u_sunlit_new) #conductance of sunlit leaves of overstorey for CO2


      Ci_u_shaded_old = Ci_u_shaded_new
      Cs_u_shaded_old = Cs_u_shaded_new
      Gs_u_shaded_old = Gs_u_shaded_new   # m/s
      Gw_u_shaded = 1.0 / (1.0 / G_u_a + 1.0 / G_u_b + 1.0 / Gs_u_shaded_new) #conductance of sunlit leaves of overstorey for water
      Gc_u_shaded = 1.0 / (1.0 / G_u_a + 1.4 / G_u_b + 1.6 / Gs_u_shaded_new) #conductance of sunlit leaves of overstorey for CO2


      # /***** Leaf temperatures module by L. He  *****/
      Leaf_Temperatures(Ta, slope, psychrometer, VPD_air, Cp_ca, Gw_o_sunlit, Gw_o_shaded, Gw_u_sunlit, Gw_u_shaded,
        Gww_o_sunlit, Gww_o_shaded, Gww_u_sunlit, Gww_u_shaded, Gh_o_sunlit, Gh_o_shaded, Gh_u_sunlit, Gh_u_shaded,
        Xcs_o[kkk], Xcl_o[kkk], Xcs_u[kkk], Xcl_u[kkk], radiation_o_sun, radiation_o_shaded, radiation_u_sun, radiation_u_shaded,
        &Tc_o_sunlit_new, &Tc_o_shaded_new, &Tc_u_sunlit_new, &Tc_u_shaded_new)


      H_o_sunlit = (Tc_o_sunlit_new - Ta) * rho_a * Cp_ca * Gh_o_sunlit
      H_o_shaded = (Tc_o_shaded_new - Ta) * rho_a * Cp_ca * Gh_o_shaded
      GH_o = H_o_sunlit * LAI_o_sunlit + H_o_shaded * LAI_o_shaded  # for next num aerodynamic conductance calculation

      if (fabs(Tc_o_sunlit_new - Tc_o_sunlit_old) < 0.02 && fabs(Tc_o_shaded_new - Tc_o_shaded_old) < 0.02 &&
          fabs(Tc_u_sunlit_new - Tc_u_sunlit_old) < 0.02 && fabs(Tc_u_shaded_new - Tc_u_shaded_old) < 0.02)

        break # break the iteration if results converge
      else
        if (num > 22)  #if the iteration does not converge

          Tc_o_sunlit_old = Ta
          Tc_o_shaded_old = Ta
          Tc_u_sunlit_old = Ta
          Tc_u_shaded_old = Ta
          break
        else
          Tc_o_sunlit_old = Tc_o_sunlit_new
          Tc_o_shaded_old = Tc_o_shaded_new
          Tc_u_sunlit_old = Tc_u_sunlit_new
          Tc_u_shaded_old = Tc_u_shaded_new
        end
      end
    end# end of while

    GPP_o_sunlit = Ac_o_sunlit * LAIo_sunlit
    GPP_o_shaded = Ac_o_shaded * LAIo_shaded
    GPP_u_sunlit = Ac_u_sunlit * LAIu_sunlit
    GPP_u_shaded = Ac_u_shaded * LAIu_shaded

    # /*****  Transpiration module by X. Luo  *****/
    transpiration(Tc_o_sunlit_new, Tc_o_shaded_new, Tc_u_sunlit_new, Tc_u_shaded_new, Ta, rh_air,
      Gw_o_sunlit, Gw_o_shaded, Gw_u_sunlit, Gw_u_shaded, LAIo_sunlit, LAIo_shaded, LAIu_sunlit, LAIu_shaded,
      &Trans_o[kkk], &Trans_u[kkk])

    # /*****  Evaporation and sublimation from canopy by X. Luo  *****/
    evaporation_canopy(Tc_o_sunlit_new, Tc_o_shaded_new, Tc_u_sunlit_new, Tc_u_shaded_new, Ta, rh_air,
      Gww_o_sunlit, Gww_o_shaded, Gww_u_sunlit, Gww_u_shaded,
      LAI_o_sunlit, LAI_o_shaded, LAI_u_sunlit, LAI_u_shaded,
      Xcl_o[kkk], Xcl_u[kkk], Xcs_o[kkk], Xcs_u[kkk],
      &Eil_o[kkk], &Eil_u[kkk], &EiS_o[kkk], &EiS_u[kkk])

    # /*****  Rainfall stage 2 by X. Luo  *****/
    rainfall_stage2(Eil_o[kkk], Eil_u[kkk], &Wcl_o[kkk], &Wcl_u[kkk])

    # /*****  Snow pack stage 2 by X. Luo  *****/
    snowpack_stage2(EiS_o[kkk], EiS_u[kkk], &Wcs_o[kkk], &Wcs_u[kkk])

    # /*****  Evaporation from soil module by X. Luo  *****/
    Gheat_g = 1 / ra_g
    mass_water_g = rho_w * Zp

    evaporation_soil(temp_grd, Ts0[kkk-1], rh_air, radiation_g, Gheat_g, &Xg_snow[kkk],
      &Zp, &Zsp, &mass_water_g, &Wg_snow[kkk], rho_snow[kkk], soilp.thetam_prev[0], soilp.fei[0],
      &Evap_soil[kkk], &Evap_SW[kkk], &Evap_SS[kkk])


    # /*****  Soil Thermal Conductivity module by L. He  *****/
    UpdateSoilThermalConductivity(soilp)
    Update_Cs(soilp)

    # /*****  Surface temperature by X. Luo  *****/

    Cs[0][kkk] = soilp.Cs[0]  # added
    Cs[1][kkk] = soilp.Cs[0]
    Tc_u[kkk] = Tcu# added
    lambda[1] = soilp.lambda[0]
    d_soil[1] = soilp.d_soil[0]
    Tm[1][kkk-1] = soilp.temp_soil_p[1] # first place is temp_soil_p[0]?
    Tm[0][kkk-1] = soilp.temp_soil_p[0]
    G[1][kkk] = soilp.G[0]

    surface_temperature(Ta, rh_air, Zsp, Zp,
      Cs[1][kkk], Cs[0][kkk], Gheat_g, d_soil[1], rho_snow[kkk], Tc_u[kkk],
      radiation_g, Evap_soil[kkk], Evap_SW[kkk], Evap_SS[kkk],
      lambda[1], Xg_snow[kkk], G[1][kkk],
      Ts0[kkk-1], Tm[1][kkk-1], Tm[0][kkk-1], Tsn0[kkk-1],
      Tsm0[kkk-1], Tsn1[kkk-1], Tsn2[kkk-1],
      &Ts0[kkk], &Tm[0][kkk], &Tsn0[kkk],
      &Tsm0[kkk], &Tsn1[kkk], &Tsn2[kkk],
      &G[0][kkk])

    soilp.temp_soil_c[0] = Tm[0][kkk]

    # /*****  Snow Pack Stage 3 module by X. Luo  *****/
    snowpack_stage3(Ta, Tsn0[kkk], Tsn0[kkk-1], rho_snow[kkk], &Zsp, &Zp, &Wg_snow[kkk])

    # /*****  Sensible heat flux module by X. Luo  *****/
    sensible_heat(Tc_o_sunlit_new, Tc_o_shaded_new, Tc_u_sunlit_new, Tc_u_shaded_new, Ts0[kkk], Ta, rh_air,
      Gh_o_sunlit, Gh_o_shaded, Gh_u_sunlit, Gh_u_shaded, Gheat_g,
      LAI_o_sunlit, LAI_o_shaded, LAI_u_sunlit, LAI_u_shaded,
      &Qhc_o[kkk], &Qhc_u[kkk], &Qhg[kkk])


    # /*****  Soil water module by L. He  *****/
    # in-process value check
    #if(jday==75 && rstep==8 && kkk==5)
    #if(jday==75 && rstep==8 )
    #printf("%d, %d, %f, %f\n", jday, rstep, soilp.thetam_prev[0], soilp.f_soilwater);

    soilp.Zsp = Zsp
    soilp.G[0] = G[0][kkk]

    UpdateHeatFlux(soilp, Xg_snow[kkk], lambda_snow[kkk], Tsn0[kkk], Ta, kstep)

    Soil_Water_Uptake(soilp, Trans_o[kkk], Trans_u[kkk], Evap_soil[kkk])

    soilp.r_rain_g = r_rain_g[kkk]
    soilp.Zp = Zp

    UpdateSoilMoisture(soilp, kstep)
    Zp = soilp.Zp

  end  # The end of kkk loop


  kkk = kloop # the last step

  Tsn1[kkk] = clamp(Tsn1[kkk], -40.0, 40.0)
  Tsn2[kkk] = clamp(Tsn2[kkk], -40.0, 40.0)
  # if (Tsn1[kkk]>40) Tsn2[kkk]=40;
  # if (Tsn1[kkk]<-40) Tsn2[kkk]=-40;

  # if (Tsn2[kkk]>40) Tsn2[kkk]=40;
  # if (Tsn2[kkk]<-40) Tsn2[kkk]=-40;

  var_n[3+1] = Ts0[kkk]        # To: The temperature of ground surface
  var_n[4+1] = Tsn0[kkk]       # To: The temperature of ground surface
  var_n[5+1] = Tsm0[kkk]
  var_n[6+1] = Tsn1[kkk]       # To: The temperature of ground surface
  var_n[7+1] = Tsn2[kkk]       # To: The temperature of ground surface
  var_n[11+1] = Qhc_o[kkk]
  var_n[15+1] = Wcl_o[kkk]
  var_n[16+1] = Wcs_o[kkk]     # the mass of intercepted liquid water and snow, overstory
  var_n[18+1] = Wcl_u[kkk]
  var_n[19+1] = Wcs_u[kkk]     # the mass of intercepted liquid water and snow, overstory
  var_n[20+1] = Wg_snow[kkk]   # the fraction of ground surface covered by snow and snow mass

  Lv_liquid = (2.501 - 0.00237 * Ta) * 1000000  # The latent heat of water vaporization in j/kg

  mid_res.Net_Rad = radiation_o + radiation_u + radiation_g
  mid_res.LH = Lv_liquid * (Trans_o[kkk] + Eil_o[kkk] + Trans_u[kkk] + Eil_u[kkk] + Evap_soil[kkk] + Evap_SW[kkk]) +
               Lv_solid * (EiS_o[kkk] + EiS_u[kkk] + Evap_SS[kkk])
  #  LH:  total latent heat flux

  mid_res.SH = Qhc_o[kkk] + Qhc_u[kkk] + Qhg[kkk] # SH: total sensible heat flux

  mid_res.Trans = (Trans_o[kkk] + Trans_u[kkk]) * step# total transpiration  mm/step
  mid_res.Evap = (Eil_o[kkk] + Eil_u[kkk] + Evap_soil[kkk] + Evap_SW[kkk] + EiS_o[kkk] + EiS_u[kkk] + Evap_SS[kkk]) * step   # total evaporation . mm/step

  mid_res.gpp_o_sunlit = GPP_o_sunlit   # umol C/m2/s
  mid_res.gpp_u_sunlit = GPP_u_sunlit
  mid_res.gpp_o_shaded = GPP_o_shaded
  mid_res.gpp_u_shaded = GPP_u_shaded

  # total GPP . gC/m2/step
  mid_res.GPP = (GPP_o_sunlit + GPP_o_shaded + GPP_u_sunlit + GPP_u_shaded) * 12 * step * 0.000001
  # return
  nothing
end
