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
  lai::T, clumping::T, param::Vector{T}, meteo::ClimateData, CosZs::T,
  var_o::Vector{T}, var_n::Vector{T}, soilp::Union{Soil_c,Soil},
  Ra::Radiation,
  mid_res::Results, mid_ET::OutputET, var::InterTempVars; kw...) where {T}

  # var = var2
  # var = InterTempVars()
  init_vars!(var)
  # reset!(var.TempLeafs)
  @unpack Cc_new, Cs_old, Cs_new, Ci_old, 
    Tc_old, Tc_new, Gs_old, Gc, Gh, Gw, Gww, 
    Gs_new, Ac, Ci_new, Rn, Rns, Rnl, 
    leleaf, GPP, LAI, PAI = var.TempLeafs

  # Cc_new = Leaf()
  # Cs_old = Leaf()
  # Cs_new = Leaf()
  # Ci_old = Leaf()
  # Tc_old = Leaf()
  # Tc_new = Leaf()
  # Gs_old = Leaf()
  # # to the reference height above the canopy
  # Gc = Leaf()  # total conductance for CO2 from the intercellular space of the leaves
  # Gh = Leaf()  # total conductance for heat transfer from the leaf surface 
  # Gw = Leaf()  # total conductance for water from the intercellular space of the leaves
  # Gww = Leaf() # total conductance for water from the surface of the leaves

  # Gs_new = Leaf()
  # Ac = Leaf()
  # Ci_new = Leaf()

  # Rn = Leaf()
  # Rns = Leaf()
  # Rnl = Leaf()
  # leleaf = Leaf()
  # GPP = Leaf()
  # LAI = Leaf()
  # PAI = Leaf()

  dz = zeros(layer + 1)
  lambda = zeros(layer + 2)

  ra_o = 0.0
  ra_u = 0.0
  ra_g = 0.0
  Ga_o = 0.0
  Gb_o = 0.0
  Ga_u = 0.0
  Gb_u = 0.0

  radiation_o = 0.0
  radiation_u = 0.0
  radiation_g = 0.0

  Tco = 0.0
  Tcu = 0.0

  alpha_sat = param[24+1]       # albedo of saturated/dry soil for module rainfall 1
  alpha_dry = param[25+1]       # the albedo of dry soil
  canopyh_o = param[29+1]       # to be used for module aerodynamic_conductance
  canopyh_u = param[30+1]
  height_wind_sp = param[31+1]  # the_height_to_measure_wind_speed, for module aerodynamic_conductance
  m_h2o = param[33+1]           # to be used for module photosynthesis
  b_h2o = param[34+1]

  # /*****  Vcmax-Nitrogen calculations，by G.Mo，Apr. 2011  *****/
  if (CosZs > 0) # day time
    # parameters for Vcmax-Nitrogen calculations
    G_theta = 0.5 # assuming a spherical leaf angle distribution
    Kn = 0.3      # 0.713/2.4
    K = G_theta * clumping / CosZs
    Vcmax0 = param[36+1]

    expr1 = 1 - exp(-K * lai)
    expr2 = 1 - exp(-lai * (Kn + K))
    expr3 = 1 - exp(-Kn * lai)

    # Formulas based on Chen et al., 2012, GBC
    if (expr1 > 0)
      Vcmax_sunlit = Vcmax0 * param[47+1] * param[46+1] * K * expr2 / (Kn + K) / expr1
    else
      Vcmax_sunlit = Vcmax0
    end

    if (K > 0 && lai > expr1 / K)
      Vcmax_shaded = Vcmax0 * param[47+1] * param[46+1] * (expr3 / Kn - expr2 / (Kn + K)) / (lai - expr1 / K)
    else
      Vcmax_shaded = Vcmax0
    end
  end

  # /*****  LAI calculation module, by B. Chen  *****/
  lai_o = lai < 0.1 ? 0.1 : lai
  
  landcover = Int(param[4+1])
  if (landcover == 25 || landcover == 40)
    lai_u = 0.01
  else
    lai_u = 1.18 * exp(-0.99 * lai_o)
  end
  (lai_u > lai_o) && (lai_u = 0.01)
  
  stem_o = param[8+1] * 0.2    # parameter[8].LAI max overstory
  stem_u = param[9+1] * 0.2    # parameter[9].LAI max understory

  # lai2: separate lai into sunlit and shaded portions
  lai2(clumping, CosZs, stem_o, stem_u, lai_o, lai_u, LAI, PAI)

  # /*****  Initialization of this time step  *****/
  Rs      = meteo.Srad
  rh_air  = meteo.rh
  wind_sp = meteo.wind
  prcp    = meteo.rain / step  # precipitation in meters
  Ta      = meteo.temp

  es = cal_es(Ta)  # to estimate saturated water vapor pressure in kpa
  ea = es * rh_air / 100                   # used in `photosynthesis`
  VPD = es - ea                            # water vapor deficit at the reference height

  q_ca = 0.622 * ea / (101.35 - 0.378 * ea)  # in g/g, unitless
  cp = Cpd * (1 + 0.84 * q_ca)
  slope = cal_slope(Ta)
  gamma = 0.066

  if (Rs <= 0)
    alpha_v_o = 0.0
    alpha_n_o = 0.0
    alpha_v_u = 0.0
    alpha_n_u = 0.0
  else
    alpha_v_o = param[22+1]
    alpha_n_o = param[23+1]
    alpha_v_u = param[22+1]
    alpha_n_u = param[23+1]
  end

  # Ground surface temperature
  var.Ts0[1] = clamp(var_o[3+1], Ta - 2.0, Ta + 2.0)  # ground0
  var.Tsm0[1] = clamp(var_o[5+1], Ta - 2.0, Ta + 2.0) # 
  var.Tsn0[1] = clamp(var_o[4+1], Ta - 2.0, Ta + 2.0) # snow0
  var.Tsn1[1] = clamp(var_o[6+1], Ta - 2.0, Ta + 2.0) # snow1
  var.Tsn2[1] = clamp(var_o[7+1], Ta - 2.0, Ta + 2.0) # snow2

  var.Qhc_o[1] = var_o[11+1]
  var.Wcl_o[1] = var_o[15+1]
  var.Wcs_o[1] = var_o[16+1]    # the mass of intercepted liquid water and snow, overstory 

  # the evaporation rate of rain and snow--in kg/m^2/s, understory
  var.Wcl_u[1] = var_o[18+1]
  var.Wcs_u[1] = var_o[19+1]    # the mass of intercepted liquid water and snow, overstory
  var.Wg_snow[1] = var_o[20+1]  # thr fraction of ground surface covered in snow and snow mass

  mass_water_g = Ref(0.0)
  Zsp = Ref(soilp.Zsp)
  Zp = Ref(soilp.Zp)
  (Zp[] < 0.001) && (Zp[] = 0.0)
  
  init_leaf_dbl(Tc_old, Ta - 0.5)
  
  # /*****  Ten time intervals in a hourly time step.6min or 360s per loop  ******/
  @inbounds for kkk = 2:kloop+1
    
    var.Xcs_o[kkk], var.Xcs_u[kkk], var.Xg_snow[kkk] = snowpack_stage1_jl(Ta, prcp,
      # var.Wcs_o[kkk-1], var.Wcs_u[kkk-1], var.Wg_snow[kkk-1],
      Ref(var.Wcs_o, kkk), Ref(var.Wcs_u, kkk), Ref(var.Wg_snow, kkk),
      lai_o, lai_u, clumping,
      Ref(var.Ac_snow_o, kkk), Ref(var.Ac_snow_u, kkk),
      Ref(var.rho_snow, kkk), Zsp,
      Ref(var.alpha_v_sw, kkk), Ref(var.alpha_n_sw, kkk)) # by X. Luo 

    # /*****  Rain fall stage 1 by X. Luo  *****/
    var.Wcl_o[kkk], var.Wcl_u[kkk], var.Xcl_o[kkk], var.Xcl_u[kkk], var.r_rain_g[kkk] = 
      rainfall_stage1_jl(Ta, prcp, var.Wcl_o[kkk-1], var.Wcl_u[kkk-1], lai_o, lai_u, clumping)

    # Old version
    # if(θ[0][kkk-1]<soilp.theta_vwp[1]*0.5) var.alpha_g = alpha_dry;
    # else var.alpha_g = (θ[0][kkk-1]-soilp.theta_vwp[1]*0.5)/(soilp.θ_sat[1]-soilp.theta_vwp[1]*0.5) * (alpha_sat - alpha_dry) + alpha_dry;
    if (soilp.θ_prev[2] < soilp.theta_vwp[2] * 0.5)
      alpha_g = alpha_dry
    else
      alpha_g = (soilp.θ_prev[2] - soilp.theta_vwp[2] * 0.5) / (soilp.θ_sat[2] - soilp.theta_vwp[2] * 0.5) * (alpha_sat - alpha_dry) + alpha_dry
    end

    alpha_v_g = 2.0 / 3.0 * alpha_g
    alpha_n_g = 4.0 / 3.0 * alpha_g

    # /*****  Soil water factor module by L. He  *****/
    soil_water_factor_v2(soilp)
    f_soilwater = min(soilp.f_soilwater, 1.0) # used in `photosynthesis`

    GH_o = var.Qhc_o[kkk-1]# to be used as the init. for module aerodynamic_conductance

    init_leaf_dbl(Ci_old, 0.7 * CO2_air)
    init_leaf_dbl2(Gs_old, 1.0 / 200.0, 1.0 / 300.0)

    percArea_snow_o = var.Ac_snow_o[kkk] / lai_o / 2
    percArea_snow_u = var.Ac_snow_u[kkk] / lai_u / 2

    temp_grd = Ta   # ground temperature substituted by air temperature

    num = 0
    while true # iteration for BWB equation until results converge
      num = num + 1

      # /***** Aerodynamic_conductance module by G.Mo  *****/
      ra_o, ra_u, ra_g, Ga_o, Gb_o, Ga_u, Gb_u =
        aerodynamic_conductance_jl(canopyh_o, canopyh_u, height_wind_sp, clumping, Ta, wind_sp, GH_o,
          lai_o + stem_o, lai_u + stem_u)

      init_leaf_dbl2(Gh,
        1.0 / (1.0 / Ga_o + 0.5 / Gb_o),
        1.0 / (1.0 / Ga_u + 0.5 / Gb_u))
      init_leaf_dbl2(Gww,
        1.0 / (1.0 / Ga_o + 1.0 / Gb_o + 100),
        1.0 / (1.0 / Ga_u + 1.0 / Gb_u + 100))

      # temperatures of overstorey and understorey canopies
      Tco = (Tc_old.o_sunlit * PAI.o_sunlit + Tc_old.o_shaded * PAI.o_shaded) / (PAI.o_sunlit + PAI.o_shaded)
      Tcu = (Tc_old.u_sunlit * PAI.u_sunlit + Tc_old.u_shaded * PAI.u_shaded) / (PAI.u_sunlit + PAI.u_shaded)

      # /*****  Net radiation at canopy and leaf level module by X.Luo  *****/
      radiation_o, radiation_u, radiation_g = netRadiation_jl(Rs, CosZs, Tco, Tcu, temp_grd,
        lai_o, lai_u, lai_o + stem_o, lai_u + stem_u, PAI,
        clumping, Ta, rh_air, var.alpha_v_sw[kkk], var.alpha_n_sw[kkk],
        percArea_snow_o, percArea_snow_u,
        var.Xg_snow[kkk],
        alpha_v_o, alpha_n_o, alpha_v_u, alpha_n_u,
        alpha_v_g, alpha_n_g, Rn, Rns, Rnl, Ra)

      # /*****  Photosynthesis module by B. Chen  *****/
      update_Gw!(Gw, Gs_old, Ga_o, Ga_u, Gb_o, Gb_u) # conductance for water
      latent_heat!(leleaf, Gw, VPD, slope, Tc_old, Ta, ρₐ, cp, gamma)

      if (CosZs > 0)
        photosynthesis(Tc_old, Rns, Ci_old, leleaf,
          Ta, ea, f_soilwater, b_h2o, m_h2o,
          Gb_o, Gb_u, Vcmax_sunlit, Vcmax_shaded,
          Gs_new, Ac, Ci_new; version="julia")
      else
        init_leaf_dbl(Gs_new, 0.0001)
        init_leaf_dbl(Ac, 0.0)
        init_leaf_dbl(Ci_new, CO2_air * 0.7)

        init_leaf_dbl(Cs_new, CO2_air)              # not used
        init_leaf_dbl(Cc_new, CO2_air * 0.7 * 0.8)  # not used
      end

      init_leaf_struct(Ci_old, Ci_new)
      init_leaf_struct(Cs_old, Cs_new)
      init_leaf_struct(Gs_old, Gs_new)

      update_Gw!(Gw, Gs_new, Ga_o, Ga_u, Gb_o, Gb_u)
      update_Gc!(Gc, Gs_new, Ga_o, Ga_u, Gb_o, Gb_u)

      # /***** Leaf temperatures module by L. He  *****/
      Leaf_Temperatures_jl(Ta, slope, gamma, VPD, cp,
        Gw, Gww, Gh,
        var.Xcs_o[kkk], var.Xcl_o[kkk], var.Xcs_u[kkk], var.Xcl_u[kkk],
        Rn, Tc_new)

      H_o_sunlit = (Tc_new.o_sunlit - Ta) * ρₐ * cp * Gh.o_sunlit
      H_o_shaded = (Tc_new.o_shaded - Ta) * ρₐ * cp * Gh.o_shaded
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

    var.Trans_o[kkk], var.Trans_u[kkk] = transpiration_jl(Tc_new, Ta, rh_air, Gw, LAI) # X. Luo

    # /*****  Evaporation and sublimation from canopy by X. Luo  *****/
    var.Eil_o[kkk], var.Eil_u[kkk], var.EiS_o[kkk], var.EiS_u[kkk] = evaporation_canopy_jl(
      Tc_new, Ta, rh_air, 
      Gww, PAI,
      var.Xcl_o[kkk], var.Xcl_u[kkk], var.Xcs_o[kkk], var.Xcs_u[kkk])

    rainfall_stage2_jl(var.Eil_o[kkk], var.Eil_u[kkk], Ref(var.Wcl_o, kkk), Ref(var.Wcl_u, kkk)) # X. Luo

    snowpack_stage2_jl(var.EiS_o[kkk], var.EiS_u[kkk], Ref(var.Wcs_o, kkk), Ref(var.Wcs_u, kkk)) # X. Luo

    # /*****  Evaporation from soil module by X. Luo  *****/
    Gheat_g = 1 / ra_g
    mass_water_g[] = ρ_w * Zp[]

    var.Evap_soil[kkk], var.Evap_SW[kkk], var.Evap_SS[kkk] = 
      evaporation_soil_jl(temp_grd, var.Ts0[kkk-1], rh_air, radiation_g, Gheat_g,
        Ref(var.Xg_snow, kkk), Zp, Zsp, mass_water_g, Ref(var.Wg_snow, kkk), # Ref
        var.rho_snow[kkk], soilp.θ_prev[1], soilp.θ_sat[1])
        # Ref(var.Evap_soil, kkk), Ref(var.Evap_SW, kkk), Ref(var.Evap_SS, kkk)

    # /*****  Soil Thermal Conductivity module by L. He  *****/
    UpdateSoilThermalConductivity(soilp)
    Update_Cs(soilp)

    # /*****  Surface temperature by X. Luo  *****/
    var.Cs .= 0.0
    var.Tm .= 0.0
    var.G .= 0.0

    var.Cs[1, kkk] = soilp.Cs[1]  # added
    var.Cs[2, kkk] = soilp.Cs[1]
    var.Tc_u[kkk]  = Tcu           # added
    lambda[2]      = soilp.lambda[1]
    dz[2]      = soilp.dz[1]

    var.Tm[1, kkk-1] = soilp.Tsoil_p[1]
    var.Tm[2, kkk-1] = soilp.Tsoil_p[2] # first place is Tsoil_p[0]?
    var.G[2, kkk] = soilp.G[1]

    var.G[1, kkk], var.Ts0[kkk], var.Tm[1, kkk], var.Tsm0[kkk],
    var.Tsn0[kkk], var.Tsn1[kkk], var.Tsn2[kkk] =
      surface_temperature_jl(Ta, rh_air, Zsp[], Zp[],
        var.Cs[2, kkk], var.Cs[1, kkk], Gheat_g, dz[2], var.rho_snow[kkk], var.Tc_u[kkk],
        radiation_g, var.Evap_soil[kkk], var.Evap_SW[kkk], var.Evap_SS[kkk],
        lambda[2], 
        var.Xg_snow[kkk], var.G[2, kkk], 
        var.Ts0[kkk-1], 
        # T_soil1_last::FT, T_any0_last::FT, T_soil0_last::FT,
        var.Tm[2, kkk-1], var.Tm[1, kkk-1], var.Tsm0[kkk-1], 
        var.Tsn0[kkk-1], var.Tsn1[kkk-1], var.Tsn2[kkk-1])

    Update_Tsoil_c(soilp, var.Tm[1, kkk])
    # soilp.Tsoil_c[1] = var.Tm[1, kkk]

    snowpack_stage3_jl(Ta, var.Tsn0[kkk], var.Tsn0[kkk-1], var.rho_snow[kkk], Zsp, Zp, Ref(var.Wg_snow, kkk)) # X. Luo

    var.Qhc_o[kkk], var.Qhc_u[kkk], var.Qhg[kkk] =
      sensible_heat_jl(Tc_new, var.Ts0[kkk], Ta, rh_air,
        Gh, Gheat_g, PAI) # X. Luo

    # /*****  Soil water module by L. He  *****/
    soilp.Zsp = Zsp[]
    Update_G(soilp, var.G[1, kkk])
    # soilp.G[1] = G[1][kkk]

    UpdateHeatFlux(soilp, var.Xg_snow[kkk], var.lambda_snow[kkk], var.Tsn0[kkk], Ta, kstep)
    Soil_Water_Uptake(soilp, var.Trans_o[kkk], var.Trans_u[kkk], var.Evap_soil[kkk])

    soilp.r_rain_g = var.r_rain_g[kkk]
    soilp.Zp = Zp[]

    UpdateSoilMoisture(soilp, kstep)
    Zp[] = soilp.Zp
  end  # The end of kkk loop

  kkk = kloop + 1# the last step
  var.Tsn1[kkk] = clamp(var.Tsn1[kkk], -40.0, 40.0)
  var.Tsn2[kkk] = clamp(var.Tsn2[kkk], -40.0, 40.0)

  var_n[3+1] = var.Ts0[kkk]        # To: The temperature of ground surface
  var_n[4+1] = var.Tsn0[kkk]       # To: The temperature of ground surface
  var_n[5+1] = var.Tsm0[kkk]
  var_n[6+1] = var.Tsn1[kkk]       # To: The temperature of ground surface
  var_n[7+1] = var.Tsn2[kkk]       # To: The temperature of ground surface
  var_n[11+1] = var.Qhc_o[kkk]
  var_n[15+1] = var.Wcl_o[kkk]
  var_n[16+1] = var.Wcs_o[kkk]     # the mass of intercepted liquid water and snow, overstory
  var_n[18+1] = var.Wcl_u[kkk]
  var_n[19+1] = var.Wcs_u[kkk]     # the mass of intercepted liquid water and snow, overstory
  var_n[20+1] = var.Wg_snow[kkk]   # the fraction of ground surface covered by snow and snow mass

  mid_res.Net_Rad = radiation_o + radiation_u + radiation_g

  OutputET!(mid_ET,
    var.Trans_o, var.Trans_u,
    var.Eil_o, var.Eil_u,
    var.EiS_o, var.EiS_u,
    var.Evap_soil, var.Evap_SW, var.Evap_SS, var.Qhc_o, var.Qhc_u, var.Qhg, kkk)
  update_ET!(mid_ET, mid_res, Ta)

  mid_res.gpp_o_sunlit = GPP.o_sunlit   # umol C/m2/s
  mid_res.gpp_u_sunlit = GPP.u_sunlit
  mid_res.gpp_o_shaded = GPP.o_shaded
  mid_res.gpp_u_shaded = GPP.u_shaded

  # total GPP . gC/m2/step
  mid_res.GPP = (GPP.o_sunlit + GPP.o_shaded + GPP.u_sunlit + GPP.u_shaded) * 12 * step * 0.000001
  nothing
end 
