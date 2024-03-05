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
  var_o::Vector{T}, var_n::Vector{T}, soil::AbstractSoil,
  Ra::Radiation,
  mid_res::Results, mid_ET::OutputET, var::InterTempVars; kw...) where {T}

  # var = InterTempVars()
  init_vars!(var)
  # reset!(var.TempLeafs)
  @unpack Cc_new, Cs_old, Cs_new, Ci_old,
  Tc_old, Tc_new, Gs_old, Gc, Gh, Gw, Gww,
  Gs_new, Ac, Ci_new, Rn, Rns, Rnl,
  leleaf, GPP, LAI, PAI = var.TempLeafs

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

  # /*****  Vcmax-Nitrogen calculations，by G.Mo，Apr. 2011  *****/
  Vcmax_sunlit, Vcmax_shaded = VCmax(lai, clumping, CosZs, param)

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
  Rs = meteo.Srad
  RH = meteo.rh
  wind = meteo.wind
  prcp = meteo.rain / step  # precipitation in meters
  Ta = meteo.temp

  met = meteo_pack_jl(Ta, RH)
  (; Δ, γ, cp, VPD, ea) = met

  α_sat = param[24+1]       # albedo of saturated/dry soil for module rainfall 1
  α_dry = param[25+1]       # the albedo of dry soil
  canopyh_o = param[29+1]       # to be used for module aerodynamic_conductance
  canopyh_u = param[30+1]
  height_wind_sp = param[31+1]  # the_height_to_measure_wind_speed, for module aerodynamic_conductance
  m_h2o = param[33+1]           # to be used for module photosynthesis
  b_h2o = param[34+1]

  if (Rs <= 0)
    α_v_o = 0.0
    α_n_o = 0.0
    α_v_u = 0.0
    α_n_u = 0.0
  else
    α_v_o = param[22+1]
    α_n_o = param[23+1]
    α_v_u = param[22+1]
    α_n_u = param[23+1]
  end

  init_leaf_dbl(Tc_old, Ta - 0.5)

  # Ground surface temperature
  var.Ts0[1] = clamp(var_o[3+1], Ta - 2.0, Ta + 2.0)  # ground0
  var.Tsm0[1] = clamp(var_o[5+1], Ta - 2.0, Ta + 2.0) # 
  var.Tsn0[1] = clamp(var_o[4+1], Ta - 2.0, Ta + 2.0) # snow0
  var.Tsn1[1] = clamp(var_o[6+1], Ta - 2.0, Ta + 2.0) # snow1
  var.Tsn2[1] = clamp(var_o[7+1], Ta - 2.0, Ta + 2.0) # snow2

  var.Qhc_o[1] = var_o[11+1]

  depth_snow = soil.depth_snow
  depth_water = soil.depth_water
  (depth_water < 0.001) && (depth_water = 0.0)

  # the evaporation rate of rain and snow--in kg/m^2/s, 
  # understory the mass of intercepted liquid water and snow, overstory
  f_snow = Layer3(0.0) # perc_snow
  m_snow = Layer3(0.0) # mass_snow
  m_snow_pre = Layer3(var_o[16+1], var_o[19+1], var_o[20+1]) # mass_snow

  # the mass of intercepted liquid water and snow, overstory 
  f_water = Layer2() # perc_water
  m_water = Layer2() # mass_water
  m_water_pre = Layer2(var_o[15+1], var_o[18+1])
  A_snow = Layer2()

  ρ_snow = init_dbl()
  α_v_sw = init_dbl()
  α_n_sw = init_dbl()

  # /*****  Ten time intervals in a hourly time step.6min or 360s per loop  ******/
  @inbounds for k = 2:kloop+1
    ρ_snow[] = 0.0 # TODO: debug C, might error
    α_v_sw[], α_n_sw[] = 0.0, 0.0

    depth_snow = snowpack_stage1_jl(Ta, prcp,
      lai_o, lai_u,
      clumping,
      m_snow_pre, m_snow, f_snow, A_snow,
      depth_snow, ρ_snow,
      α_v_sw, α_n_sw) # by X. Luo 

    # /*****  Rain fall stage 1 by X. Luo  *****/
    var.r_rain_g[k] = rainfall_stage1_jl(Ta, prcp, f_water, m_water, m_water_pre, lai_o, lai_u, clumping)

    if (soil.θ_prev[2] < soil.θ_vwp[2] * 0.5)
      α_g = α_dry
    else
      α_g = (soil.θ_prev[2] - soil.θ_vwp[2] * 0.5) / (soil.θ_sat[2] - soil.θ_vwp[2] * 0.5) * (α_sat - α_dry) + α_dry
    end

    α_v_g = 2.0 / 3.0 * α_g
    α_n_g = 4.0 / 3.0 * α_g

    # /*****  Soil water factor module by L. He  *****/
    soil_water_factor_v2(soil)
    f_soilwater = min(soil.f_soilwater, 1.0) # used in `photosynthesis`

    GH_o = var.Qhc_o[k-1]# to be used as the init. for module aerodynamic_conductance

    init_leaf_dbl(Ci_old, 0.7 * CO2_air)
    init_leaf_dbl2(Gs_old, 1.0 / 200.0, 1.0 / 300.0)

    percArea_snow_o = A_snow.o / lai_o / 2
    percArea_snow_u = A_snow.u / lai_u / 2

    temp_grd = Ta   # ground temperature substituted by air temperature

    num = 0
    while true # iteration for BWB equation until results converge
      num = num + 1
      # /***** Aerodynamic_conductance module by G.Mo  *****/
      ra_o, ra_u, ra_g, Ga_o, Gb_o, Ga_u, Gb_u =
        aerodynamic_conductance_jl(canopyh_o, canopyh_u, height_wind_sp, clumping, Ta, wind, GH_o,
          lai_o + stem_o, lai_u + stem_u)

      init_leaf_dbl2(Gh,
        1.0 / (1.0 / Ga_o + 0.5 / Gb_o),
        1.0 / (1.0 / Ga_u + 0.5 / Gb_u))
      init_leaf_dbl2(Gww,
        1.0 / (1.0 / Ga_o + 1.0 / Gb_o + 100),
        1.0 / (1.0 / Ga_u + 1.0 / Gb_u + 100))

      # temperatures of overstorey and understorey canopies
      Tco = (Tc_old.o_sunlit * PAI.o_sunlit + Tc_old.o_shaded * PAI.o_shaded) /
            (PAI.o_sunlit + PAI.o_shaded)
      Tcu = (Tc_old.u_sunlit * PAI.u_sunlit + Tc_old.u_shaded * PAI.u_shaded) /
            (PAI.u_sunlit + PAI.u_shaded)

      # /*****  Net radiation at canopy and leaf level module by X.Luo  *****/
      radiation_o, radiation_u, radiation_g = netRadiation_jl(Rs, CosZs, Tco, Tcu, temp_grd,
        lai_o, lai_u, lai_o + stem_o, lai_u + stem_u, PAI,
        clumping, Ta, RH,
        α_v_sw[], α_n_sw[],
        percArea_snow_o, percArea_snow_u,
        f_snow.g,
        α_v_o, α_n_o,
        α_v_u, α_n_u,
        α_v_g, α_n_g,
        Rn, Rns, Rnl, Ra)

      # /*****  Photosynthesis module by B. Chen  *****/
      update_Gw!(Gw, Gs_old, Ga_o, Ga_u, Gb_o, Gb_u) # conductance for water
      latent_heat!(leleaf, Gw, VPD, Δ, Tc_old, Ta, ρₐ, cp, γ)

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

      set!(Ci_old, Ci_new)
      set!(Cs_old, Cs_new)
      set!(Gs_old, Gs_new)

      update_Gw!(Gw, Gs_new, Ga_o, Ga_u, Gb_o, Gb_u)
      update_Gc!(Gc, Gs_new, Ga_o, Ga_u, Gb_o, Gb_u)

      # /***** Leaf temperatures module by L. He  *****/
      Leaf_Temperatures_jl(Ta, Δ, γ, VPD, cp,
        Gw, Gww, Gh, f_water, f_snow, Rn, Tc_new)

      H_o_sunlit = (Tc_new.o_sunlit - Ta) * ρₐ * cp * Gh.o_sunlit
      H_o_shaded = (Tc_new.o_shaded - Ta) * ρₐ * cp * Gh.o_shaded
      # for next num aerodynamic conductance calculation
      GH_o = H_o_sunlit * PAI.o_sunlit + H_o_shaded * PAI.o_shaded

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
          set!(Tc_old, Tc_new)
        end
      end
    end# end of while
    multiply!(GPP, Ac, LAI)

    var.Trans_o[k], var.Trans_u[k] = transpiration_jl(Tc_new, Ta, RH, Gw, LAI) # X. Luo

    # /*****  Evaporation and sublimation from canopy by X. Luo  *****/
    var.Eil_o[k], var.Eil_u[k], var.EiS_o[k], var.EiS_u[k] = evaporation_canopy_jl(
      Tc_new, Ta, RH,
      Gww, PAI, f_water, f_snow)

    rainfall_stage2_jl(var.Eil_o[k], var.Eil_u[k], m_water) # X. Luo
    set!(m_water_pre, m_water)

    snowpack_stage2_jl(var.EiS_o[k], var.EiS_u[k], m_snow) # X. Luo

    # /*****  Evaporation from soil module by X. Luo  *****/
    Gheat_g = 1 / ra_g
    mass_water_g = ρ_w * depth_water

    var.Evap_soil[k], var.Evap_SW[k], var.Evap_SS[k], depth_water, depth_snow =
      evaporation_soil_jl(temp_grd, var.Ts0[k-1], RH, radiation_g, Gheat_g,
        f_snow,
        depth_water, depth_snow, mass_water_g, m_snow, # Ref
        ρ_snow[], soil.θ_prev[1], soil.θ_sat[1])

    # /*****  Soil Thermal Conductivity module by L. He  *****/
    UpdateSoilThermalConductivity(soil)
    Update_Cs(soil)

    # /*****  Surface temperature by X. Luo  *****/
    var.Cs .= 0.0
    var.Tm .= 0.0
    var.G .= 0.0

    var.Cs[1, k] = soil.Cs[1]
    var.Cs[2, k] = soil.Cs[1]
    var.Tc_u[k] = Tcu
    lambda[2] = soil.lambda[1]
    dz[2] = soil.dz[1]

    var.Tm[1, k-1] = soil.Tsoil_p[1]
    var.Tm[2, k-1] = soil.Tsoil_p[2] # first place is Tsoil_p[0]?
    var.G[2, k] = soil.G[1]

    var.G[1, k], var.Ts0[k], var.Tm[1, k], var.Tsm0[k],
    var.Tsn0[k], var.Tsn1[k], var.Tsn2[k] =
      surface_temperature_jl(Ta, RH, depth_snow, depth_water,
        var.Cs[2, k], var.Cs[1, k], Gheat_g, dz[2], ρ_snow[], var.Tc_u[k],
        radiation_g, var.Evap_soil[k], var.Evap_SW[k], var.Evap_SS[k],
        lambda[2],
        f_snow.g, var.G[2, k],
        var.Ts0[k-1],
        # T_soil1_last::FT, T_any0_last::FT, T_soil0_last::FT,
        var.Tm[2, k-1], var.Tm[1, k-1], var.Tsm0[k-1],
        var.Tsn0[k-1], var.Tsn1[k-1], var.Tsn2[k-1])

    # Update_Tsoil_c(soilp, var.Tm[1, k])
    soil.Tsoil_c[1] = var.Tm[1, k]

    depth_snow, depth_water = snowpack_stage3_jl(Ta, var.Tsn0[k], var.Tsn0[k-1],
      ρ_snow[], depth_snow, depth_water, m_snow) # X. Luo
    set!(m_snow_pre, m_snow) # update snow_pre

    var.Qhc_o[k], var.Qhc_u[k], var.Qhg[k] =
      sensible_heat_jl(Tc_new, var.Ts0[k], Ta, RH,
        Gh, Gheat_g, PAI) # X. Luo

    # /*****  Soil water module by L. He  *****/
    soil.depth_snow = depth_snow
    soil.G[1] = var.G[1, k]
    # Update_G(soilp, var.G[1, k])

    UpdateHeatFlux(soil, Ta, kstep) # f_snow.g, var.lambda_snow[k], var.Tsn0[k],
    Soil_Water_Uptake(soil, var.Trans_o[k], var.Trans_u[k], var.Evap_soil[k])

    soil.r_rain_g = var.r_rain_g[k]
    soil.depth_water = depth_water

    UpdateSoilMoisture(soil, kstep)
    depth_water = soil.depth_water
  end  # The end of k loop

  k = kloop + 1# the last step
  var.Tsn1[k] = clamp(var.Tsn1[k], -40.0, 40.0)
  var.Tsn2[k] = clamp(var.Tsn2[k], -40.0, 40.0)

  var_n[3+1] = var.Ts0[k]        # To: The temperature of ground surface
  var_n[4+1] = var.Tsn0[k]       # To: The temperature of ground surface
  var_n[5+1] = var.Tsm0[k]
  var_n[6+1] = var.Tsn1[k]       # To: The temperature of ground surface
  var_n[7+1] = var.Tsn2[k]       # To: The temperature of ground surface
  var_n[11+1] = var.Qhc_o[k]

  var_n[15+1] = m_water.o
  var_n[18+1] = m_water.u

  var_n[16+1] = m_snow.o
  var_n[19+1] = m_snow.u
  var_n[20+1] = m_snow.g

  mid_res.Net_Rad = radiation_o + radiation_u + radiation_g

  OutputET!(mid_ET,
    var.Trans_o, var.Trans_u,
    var.Eil_o, var.Eil_u,
    var.EiS_o, var.EiS_u,
    var.Evap_soil, var.Evap_SW, var.Evap_SS, var.Qhc_o, var.Qhc_u, var.Qhg, k)
  update_ET!(mid_ET, mid_res, Ta)

  mid_res.gpp_o_sunlit = GPP.o_sunlit   # umol C/m2/s
  mid_res.gpp_u_sunlit = GPP.u_sunlit
  mid_res.gpp_o_shaded = GPP.o_shaded
  mid_res.gpp_u_shaded = GPP.u_shaded
  # mid_res.z_water = depth_water
  # mid_res.z_snow = depth_snow
  # mid_res.ρ_snow = ρ_snow

  # total GPP . gC/m2/step
  mid_res.GPP = (GPP.o_sunlit + GPP.o_shaded + GPP.u_sunlit + GPP.u_shaded) * 12 * step * 0.000001
  nothing
end

# TODO: 
# - var_n: 可以重写为结构体
