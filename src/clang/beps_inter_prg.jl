function inter_prg(jday, rstep, 
  lai::T, clumping::T, parameter::Vector{T}, meteo::ClimateData, CosZs::T, 
  var_o::Vector{T}, var_n::Vector{T}, soilp::Soil, mid_res::Results) where {T<:Real}

  ccall((:inter_prg, libbeps), Cvoid,
    (Cint, Cint, Cdouble, Cdouble, Ptr{Cdouble},
      Ptr{ClimateData}, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Soil}, Ptr{Results}),
    jday, rstep, lai, clumping, parameter,
    Ref(meteo), CosZs, var_o, var_n, Ref(soilp), Ref(mid_res))
end


"""
The inter-module function between main program and modules

# Arguments

- `jday`      : day of year
- `rstep`     : hour of day
- `lai`       : leaf area index
- `clumping`  : clumping index
- `parameter` : parameter array according to land cover types
- `meteo`     : meteorological data
- `CosZs`     : cosine of solar zenith angle
- `var_o`     : temporary variables array of last time step
- `var_n`     : temporary variables array of this time step
- `soilp`     : soil coefficients according to land cover types and soil textures
- `mid_res`   : results struct
"""
function inter_prg_jl(
  jday::Int, rstep::Int, 
  lai::T, clumping::T, parameter::Vector{T}, meteo::ClimateData, CosZs::T,
  var_o::Vector{T}, var_n::Vector{T}, soilp::Soil, mid_res::Results, mid_ET::OutputET) where {T}
  
  d_soil = zeros(layer + 1)
  lambda = zeros(layer + 2)

  Cc_new    = Leaf()
  Cs_old    = Leaf()
  Cs_new    = Leaf()
  Ci_old    = Leaf()
  Tc_old    = Leaf()
  Tc_new    = Leaf()
  Gs_old    = Leaf()
  Gc        = Leaf() # total conductance for CO2 from the intercellular space of the leaves to the reference height above the canopy
  Gh        = Leaf() # total conductance for heat transfer from the leaf surface to the reference height above the canopy
  Gw        = Leaf() # total conductance for water from the intercellular space of the leaves to the reference height above the canopy
  Gww       = Leaf() # total conductance for water from the surface of the leaves to the reference height above the canopy
  
  Gs_new    = LeafRef()
  Ac        = LeafRef()
  Ci_new    = LeafRef()

  radiation = Leaf()
  R         = Leaf()
  leleaf    = Leaf()
  GPP       = Leaf()

  ra_o = init_dbl()
  ra_u = init_dbl()
  ra_g = init_dbl()
  Ga_o = init_dbl()
  Gb_o = init_dbl()
  Ga_u = init_dbl()
  Gb_u = init_dbl()

  radiation_o = init_dbl()
  radiation_u = init_dbl()
  radiation_g = init_dbl()

  Tco = 0.0
  Tcu = 0.0
  Tc_u = zeros(MAX_Loop)
  Ts0 = zeros(MAX_Loop)
  Tsn0 = zeros(MAX_Loop)
  Tsm0 = zeros(MAX_Loop)
  Tsn1 = zeros(MAX_Loop)
  Tsn2 = zeros(MAX_Loop)

  Qhc_o = zeros(MAX_Loop)
  Qhc_u = zeros(MAX_Loop)
  Qhg = zeros(MAX_Loop)

  Wcl_o = zeros(MAX_Loop)
  Wcs_o = zeros(MAX_Loop) # the masses of liquid and snow on the canopy

  Xcl_o = zeros(MAX_Loop)
  Xcs_o = zeros(MAX_Loop) # the fraction of rain and snow on the canopy

  Wcl_u = zeros(MAX_Loop)
  Wcs_u = zeros(MAX_Loop) # the masses of rain and snow on the canopy
  Wg_snow = zeros(MAX_Loop)

  Xcl_u = zeros(MAX_Loop)
  Xcs_u = zeros(MAX_Loop) # the fraction of rain and snow on the canopy

  Ac_snow_o = zeros(MAX_Loop)
  Ac_snow_u = zeros(MAX_Loop)

  Xg_snow = zeros(MAX_Loop)
  rho_snow = zeros(MAX_Loop)

  alpha_v_sw = zeros(MAX_Loop)
  alpha_n_sw = zeros(MAX_Loop)

  r_rain_g = zeros(MAX_Loop)

  Trans_o = zeros(MAX_Loop)
  Trans_u = zeros(MAX_Loop)
  Eil_o = zeros(MAX_Loop)
  EiS_o = zeros(MAX_Loop)
  Eil_u = zeros(MAX_Loop)
  EiS_u = zeros(MAX_Loop)

  lambda_snow = zeros(MAX_Loop)

  Evap_soil = zeros(MAX_Loop)
  Evap_SW = zeros(MAX_Loop)
  Evap_SS = zeros(MAX_Loop)

  # # parameters for Vcmax-Nitrogen calculations
  Kn = 0.3  #0.713/2.4
  G_theta = 0.5
  psychrometer = 0.066
  
  # double K,Vcmax0,Vcmax_sunlit,Vcmax_shaded,expr1,expr2,expr3;
  # double slope_Vcmax_N, leaf_N,Jmax_sunlit,Jmax_shaded;

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

  # lai2: separate lai into sunlit and shaded portions
  LAI = Leaf()
  PAI = Leaf()
  lai2(clumping, CosZs, stem_o, stem_u, lai_o, lai_u, LAI, PAI)

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

  Qhc_o[1] = var_o[11+1]
  Wcl_o[1] = var_o[15+1]
  Wcs_o[1] = var_o[16+1]   # /* the mass of intercepted liquid water and snow, overstory */

  # the evaporation rate of rain and snow--in kg/m^2/s, understory
  Wcl_u[1] = var_o[18+1]
  Wcs_u[1] = var_o[19+1]   # /* the mass of intercepted liquid water and snow, overstory */
  Wg_snow[1] = var_o[20+1]  # /* thr fraction of ground surface covered in snow and snow mass */

  Zsp = Ref(soilp.Zsp)
  Zp = Ref(soilp.Zp)

  if (Zp[] < 0.001)
    Zp[] = 0
  end

  init_leaf_dbl(Tc_old, Ta - 0.5);
  
  # /*****  Vcmax Jmax module by L. He  *****/
  # /*****  Ten time intervals in a hourly time step.6min or 360s per loop  ******/
  # for(kkk = 1;kkk <= kloop;kkk++)
  for kkk = 2:kloop+1
    # /*****  Snow pack stage 1 by X. Luo  *****/
    snowpack_stage1(Ta, precip, 
      Wcs_o[kkk-1], Wcs_u[kkk-1], Wg_snow[kkk-1],
      Ref(Wcs_o, kkk), Ref(Wcs_u, kkk), Ref(Wg_snow, kkk), 
      lai_o, lai_u, clumping,
      Ref(Ac_snow_o, kkk), Ref(Ac_snow_u, kkk), 
      Ref(Xcs_o, kkk), Ref(Xcs_u, kkk), Ref(Xg_snow, kkk),
      Ref(rho_snow, kkk), Zsp,
      Ref(alpha_v_sw, kkk), Ref(alpha_n_sw, kkk))

    # /*****  Rain fall stage 1 by X. Luo  *****/
    rainfall_stage1(Ta, precip, Wcl_o[kkk-1], Wcl_u[kkk-1],
      lai_o, lai_u, clumping, 
      Ref(Wcl_o, kkk), Ref(Wcl_u, kkk), Ref(Xcl_o, kkk),
      Ref(Xcl_u, kkk), Ref(r_rain_g, kkk))

    # Old version
    # if(thetam[0][kkk-1]<soilp.theta_vwp[1]*0.5) alpha_g = alpha_dry;
    # else alpha_g = (thetam[0][kkk-1]-soilp.theta_vwp[1]*0.5)/(soilp.fei[1]-soilp.theta_vwp[1]*0.5) * (alpha_sat - alpha_dry) + alpha_dry;
    if (soilp.thetam_prev[2] < soilp.theta_vwp[2] * 0.5)
      alpha_g = alpha_dry
    else
      alpha_g = (soilp.thetam_prev[2] - soilp.theta_vwp[2] * 0.5) / (soilp.fei[2] - soilp.theta_vwp[2] * 0.5) * (alpha_sat - alpha_dry) + alpha_dry
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

    init_leaf_dbl(Ci_old, 0.7 * CO2_air)
    init_leaf_dbl2(Gs_old, 1.0 / 200.0, 1.0 / 300.0)

    percentArea_snow_o = Ac_snow_o[kkk] / lai_o / 2
    percentArea_snow_u = Ac_snow_u[kkk] / lai_u / 2

    temp_grd = Ta   # ground temperature substituted by air temperature
    
    num = 0
    while true # iteration for BWB equation until results converge
      num = num + 1

      # /***** Aerodynamic_conductance module by G.Mo  *****/
      aerodynamic_conductance(canopyh_o, canopyh_u, height_wind_sp, clumping, Ta, wind_sp, GH_o,
        lai_o + stem_o, lai_u + stem_u, 
        ra_o, ra_u, ra_g, Ga_o, Gb_o, Ga_u, Gb_u)

      init_leaf_dbl2(Gh, 
        1.0 / (1.0 / Ga_o[] + 0.5 / Gb_o[]), 
        1.0 / (1.0 / Ga_u[] + 0.5 / Gb_u[]))
      init_leaf_dbl2(Gww,
        1.0 / (1.0 / Ga_o[] + 1.0 / Gb_o[] + 100),
        1.0 / (1.0 / Ga_u[] + 1.0 / Gb_u[] + 100))

      # temperatures of overstorey and understorey canopies
      Tco = (Tc_old.o_sunlit * PAI.o_sunlit + Tc_old.o_shaded * PAI.o_shaded) / (PAI.o_sunlit + PAI.o_shaded)
      Tcu = (Tc_old.u_sunlit * PAI.u_sunlit + Tc_old.u_shaded * PAI.u_shaded) / (PAI.u_sunlit + PAI.u_shaded)

      # /*****  Net radiation at canopy and leaf level module by X.Luo  *****/
      netRadiation(Ks, CosZs, Tco, Tcu, temp_grd, 
        lai_o, lai_u, lai_o + stem_o, lai_u + stem_u, PAI,
        clumping, Ta, rh_air, alpha_v_sw[kkk], alpha_n_sw[kkk],
        percentArea_snow_o, percentArea_snow_u,
        Xg_snow[kkk], alpha_v_o, alpha_n_o, alpha_v_u, alpha_n_u, alpha_v_g, alpha_n_g,
        radiation_o, radiation_u, radiation_g,
        radiation, R)

      # /*****  Photosynthesis module by B. Chen  *****/
      # conductance for water
      update_Gw!(Gw, Gs_old, Ga_o[], Ga_u[], Gb_o[], Gb_u[])
      latent_heat!(leleaf, Gw, VPD_air, slope, Tc_old, Ta, rho_a, Cp_ca, psychrometer)

      if (CosZs > 0)
        photosynthesis(Tc_old, R, Ci_old, leleaf, 
          Ta, e_a10, f_soilwater, b_h2o, m_h2o, 
          Gb_o[], Gb_u[], Vcmax_sunlit, Vcmax_shaded, 
          Gs_new, Ac, Ci_new)
      else
        init_leaf_dbl(Gs_new, 0.0001);
        init_leaf_dbl(Ac, 0.0);
        init_leaf_dbl(Ci_new, CO2_air * 0.7)

        init_leaf_dbl(Cs_new, CO2_air)               # not used
        init_leaf_dbl(Cc_new, CO2_air * 0.7 * 0.8);  # not used
      end

      init_leaf_struct(Ci_old, Ci_new)
      init_leaf_struct(Cs_old, Cs_new)
      init_leaf_struct(Gs_old, Gs_new)

      update_Gw!(Gw, Gs_new, Ga_o[], Ga_u[], Gb_o[], Gb_u[])
      update_Gc!(Gc, Gs_new, Ga_o[], Ga_u[], Gb_o[], Gb_u[])

      # /***** Leaf temperatures module by L. He  *****/
      Leaf_Temperatures(Ta, slope, psychrometer, VPD_air, Cp_ca,
        Gw, Gww, Gh,
        Xcs_o[kkk], Xcl_o[kkk], Xcs_u[kkk], Xcl_u[kkk],
        radiation, Tc_new)

      H_o_sunlit = (Tc_new.o_sunlit - Ta) * rho_a * Cp_ca * Gh.o_sunlit
      H_o_shaded = (Tc_new.o_shaded - Ta) * rho_a * Cp_ca * Gh.o_shaded
      GH_o = H_o_sunlit * PAI.o_sunlit + H_o_shaded * PAI.o_shaded  # for next num aerodynamic conductance calculation

      if (abs(Tc_new.o_sunlit - Tc_old.o_sunlit) < 0.02 && 
          abs(Tc_new.o_shaded - Tc_old.o_shaded) < 0.02 &&
          abs(Tc_new.u_sunlit - Tc_old.u_sunlit) < 0.02 && 
          abs(Tc_new.u_shaded - Tc_old.u_shaded) < 0.02)
        break # break the iteration if results converge
      else
        if (num > 22)  #if the iteration does not converge
          init_leaf_dbl(Tc_old, Ta)
          break
        else
          init_leaf_struct(Tc_old, Tc_new)
        end
      end
    end# end of while
    multiply!(GPP, Ac, LAI)
    
    # /*****  Transpiration module by X. Luo  *****/
    transpiration(
      Tc_new, Ta, rh_air,
      Gw, LAI,
      Ref(Trans_o, kkk), Ref(Trans_u, kkk))

    # /*****  Evaporation and sublimation from canopy by X. Luo  *****/
    evaporation_canopy(
      Tc_new, Ta, rh_air,
      Gww, PAI,
      Xcl_o[kkk], Xcl_u[kkk], Xcs_o[kkk], Xcs_u[kkk],
      Ref(Eil_o, kkk), Ref(Eil_u, kkk), Ref(EiS_o, kkk), Ref(EiS_u, kkk))

    # /*****  Rainfall stage 2 by X. Luo  *****/
    rainfall_stage2(Eil_o[kkk], Eil_u[kkk], Ref(Wcl_o, kkk), Ref(Wcl_u, kkk))

    # /*****  Snow pack stage 2 by X. Luo  *****/
    snowpack_stage2(EiS_o[kkk], EiS_u[kkk], Ref(Wcs_o, kkk), Ref(Wcs_u, kkk))

    # /*****  Evaporation from soil module by X. Luo  *****/
    Gheat_g = 1 / ra_g[]
    mass_water_g = Ref(rho_w * Zp[])
    
    evaporation_soil(temp_grd, Ts0[kkk-1], rh_air, radiation_g[], Gheat_g, 
      Ref(Xg_snow, kkk), Zp, Zsp, mass_water_g, Ref(Wg_snow, kkk), # Ref
      rho_snow[kkk], soilp.thetam_prev[1], soilp.fei[1],
      Ref(Evap_soil, kkk), Ref(Evap_SW, kkk), Ref(Evap_SS, kkk))

    # /*****  Soil Thermal Conductivity module by L. He  *****/
    UpdateSoilThermalConductivity(soilp)
    Update_Cs(soilp)

    # /*****  Surface temperature by X. Luo  *****/
    Cs = zeros(layer + 2, MAX_Loop)
    Tm = zeros(layer + 2, MAX_Loop)
    G = zeros(layer + 2, MAX_Loop)
    
    Cs[1, kkk] = soilp.Cs[1]  # added
    Cs[2, kkk] = soilp.Cs[1]
    Tc_u[kkk] = Tcu# added
    lambda[2] = soilp.lambda[1]
    d_soil[2] = soilp.d_soil[1]
    
    Tm[2, kkk-1] = soilp.temp_soil_p[2] # first place is temp_soil_p[0]?
    Tm[1, kkk-1] = soilp.temp_soil_p[1]
    G[2, kkk] = soilp.G[1]

    ## 二维数组`Ref`如何处理？
    Ts0[kkk], Tm[1, kkk], Tsn0[kkk], Tsm0[kkk], Tsn1[kkk], Tsn2[kkk], G[1, kkk] = 
      surface_temperature(Ta, rh_air, Zsp[], Zp[],
        Cs[2, kkk], Cs[1, kkk], Gheat_g, d_soil[2], rho_snow[kkk], Tc_u[kkk],
        radiation_g[], Evap_soil[kkk], Evap_SW[kkk], Evap_SS[kkk],
        lambda[2], Xg_snow[kkk], G[2, kkk],
        Ts0[kkk-1], Tm[2, kkk-1], Tm[1, kkk-1], Tsn0[kkk-1],
        Tsm0[kkk-1], Tsn1[kkk-1], Tsn2[kkk-1])

    # @show soilp.temp_soil_c, Tm[1, kkk]
    Update_temp_soil_c(soilp, Tm[1, kkk])
    # soilp.temp_soil_c[1] = Tm[1, kkk]

    # /*****  Snow Pack Stage 3 module by X. Luo  *****/
    snowpack_stage3(Ta, Tsn0[kkk], Tsn0[kkk-1], rho_snow[kkk], Zsp, Zp, Ref(Wg_snow, kkk))

    # /*****  Sensible heat flux module by X. Luo  *****/
    sensible_heat(Tc_new, Ts0[kkk], Ta, rh_air,
      Gh, Gheat_g, PAI, 
      Ref(Qhc_o, kkk), Ref(Qhc_u, kkk), Ref(Qhg, kkk))

    # /*****  Soil water module by L. He  *****/
    soilp.Zsp = Zsp[]
    Update_G(soilp, G[1, kkk])
    # soilp.G[1] = G[1][kkk]

    UpdateHeatFlux(soilp, Xg_snow[kkk], lambda_snow[kkk], Tsn0[kkk], Ta, kstep)
    Soil_Water_Uptake(soilp, Trans_o[kkk], Trans_u[kkk], Evap_soil[kkk])

    soilp.r_rain_g = r_rain_g[kkk]
    soilp.Zp = Zp[]

    UpdateSoilMoisture(soilp, kstep)
    Zp[] = soilp.Zp
  end  # The end of kkk loop

  kkk = kloop # the last step
  Tsn1[kkk] = clamp(Tsn1[kkk], -40.0, 40.0)
  Tsn2[kkk] = clamp(Tsn2[kkk], -40.0, 40.0)
  
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

  mid_res.Net_Rad = radiation_o[] + radiation_u[] + radiation_g[]
  
  OutputET!(mid_ET, 
    Trans_o, Trans_u, 
    Eil_o, Eil_u,
    EiS_o, EiS_u, 
    Evap_soil, Evap_SW, Evap_SS, Qhc_o, Qhc_u, Qhg, kkk)
  update_ET!(mid_ET, mid_res, Ta)
  
  mid_res.gpp_o_sunlit = GPP.o_sunlit   # umol C/m2/s
  mid_res.gpp_u_sunlit = GPP.u_sunlit
  mid_res.gpp_o_shaded = GPP.o_shaded
  mid_res.gpp_u_shaded = GPP.u_shaded
  
  # total GPP . gC/m2/step
  mid_res.GPP = (GPP.o_sunlit + GPP.o_shaded + GPP.u_sunlit + GPP.u_shaded) * 12 * step * 0.000001
  # return
  nothing
end
