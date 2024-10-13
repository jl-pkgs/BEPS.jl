"""
The inter-module function between main program and modules

# Arguments

- `jday`      : day of year
- `hour`     : hour of day
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
  jday::Int, hour::Int,
  lai::T, Ω::T, param::Vector{T}, meteo::Met, CosZs::T,
  state::State{T}, soil::AbstractSoil,
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
  Vcmax_sunlit, Vcmax_shaded = VCmax(lai, Ω, CosZs, param)

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
  lai2(Ω, CosZs, stem_o, stem_u, lai_o, lai_u, LAI, PAI)

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
  z_wind = param[31+1]  # the_height_to_measure_wind_speed, for module aerodynamic_conductance
  g1_w = param[33+1]           # to be used for module photosynthesis
  g0_w = param[34+1]

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
  # [Ts0, Tsn, Tsm0, Tsn1, Tsn2]
  var.Ts0[1] = clamp(state.Ts[1], Ta - 2.0, Ta + 2.0)  # ground0
  var.Tsn0[1] = clamp(state.Ts[2], Ta - 2.0, Ta + 2.0) # snow0
  var.Tsm0[1] = clamp(state.Ts[3], Ta - 2.0, Ta + 2.0) # any
  var.Tsn1[1] = clamp(state.Ts[4], Ta - 2.0, Ta + 2.0) # snow1
  var.Tsn2[1] = clamp(state.Ts[5], Ta - 2.0, Ta + 2.0) # snow2
  var.Qhc_o[1] = state.Qhc_o

  m_snow_pre = state.m_snow# mass_snow
  m_water_pre = state.m_water
  
  z_snow = soil.z_snow
  z_water = soil.z_water
  (z_water < 0.001) && (z_water = 0.0)

  # the evaporation rate of rain and snow--in kg/m^2/s, 
  # understory the mass of intercepted liquid water and snow, overstory
  f_snow = Layer3(0.0) # perc_snow
  m_snow = Layer3(0.0) # mass_snow
  A_snow = Layer2()
  
  # the mass of intercepted liquid water and snow, overstory 
  f_water = Layer2() # perc_water
  m_water = Layer2() # mass_water

  ρ_snow = init_dbl()
  α_v_sw = init_dbl()
  α_n_sw = init_dbl()

  # /*****  Ten time intervals in a hourly time step.6min or 360s per loop  ******/
  @inbounds for k = 2:kloop+1
    ρ_snow[] = 0.0 # TODO: debug C, might error
    α_v_sw[], α_n_sw[] = 0.0, 0.0

    # set!(m_snow_pre, 0.0);
    z_snow = snowpack_stage1_jl(Ta, prcp,
      lai_o, lai_u, Ω,
      m_snow_pre, m_snow, f_snow, A_snow,
      z_snow, ρ_snow,
      α_v_sw, α_n_sw) # by X. Luo 

    # /*****  Rain fall stage 1 by X. Luo  *****/
    var.r_rain_g[k] = rainfall_stage1_jl(Ta, prcp, f_water, m_water, m_water_pre, lai_o, lai_u, Ω)

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

    perc_snow_o = A_snow.o / lai_o / 2 # area
    perc_snow_u = A_snow.u / lai_u / 2 # area

    temp_grd = Ta   # ground temperature substituted by air temperature

    num = 0
    while true # iteration for BWB equation until results converge
      num = num + 1
      # /***** Aerodynamic_conductance module by G.Mo  *****/
      ra_o, ra_u, ra_g, Ga_o, Gb_o, Ga_u, Gb_u =
        aerodynamic_conductance_jl(canopyh_o, canopyh_u, z_wind, Ω, Ta, wind, GH_o,
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
        Ω, Ta, RH,
        α_v_sw[], α_n_sw[],
        perc_snow_o, perc_snow_u,
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
          Ta, ea, f_soilwater, g0_w, g1_w,
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
    mass_water_g = ρ_w * z_water

    var.Evap_soil[k], var.Evap_SW[k], var.Evap_SS[k], z_water, z_snow =
      evaporation_soil_jl(temp_grd, var.Ts0[k-1], RH, radiation_g, Gheat_g,
        f_snow,
        z_water, z_snow, mass_water_g, m_snow, # Ref
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
      surface_temperature_jl(Ta, RH, z_snow, z_water,
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

    z_snow, z_water = snowpack_stage3_jl(Ta, var.Tsn0[k], var.Tsn0[k-1],
      ρ_snow[], z_snow, z_water, m_snow) # X. Luo
    set!(m_snow_pre, m_snow) # update snow_pre

    var.Qhc_o[k], var.Qhc_u[k], var.Qhg[k] =
      sensible_heat_jl(Tc_new, var.Ts0[k], Ta, RH,
        Gh, Gheat_g, PAI) # X. Luo

    # /*****  Soil water module by L. He  *****/
    soil.z_snow = z_snow
    soil.G[1] = var.G[1, k]
    # Update_G(soilp, var.G[1, k])

    UpdateHeatFlux(soil, Ta, kstep) # f_snow.g, var.lambda_snow[k], var.Tsn0[k],
    Soil_Water_Uptake(soil, var.Trans_o[k], var.Trans_u[k], var.Evap_soil[k])

    soil.r_rain_g = var.r_rain_g[k]
    soil.z_water = z_water

    UpdateSoilMoisture(soil, kstep)
    z_water = soil.z_water
  end  # The end of k loop

  k = kloop + 1# the last step
  var.Tsn1[k] = clamp(var.Tsn1[k], -40.0, 40.0)
  var.Tsn2[k] = clamp(var.Tsn2[k], -40.0, 40.0)

  state.Ts[1] = var.Ts0[k]
  state.Ts[2] = var.Tsn0[k]
  state.Ts[3] = var.Tsm0[k]
  state.Ts[4] = var.Tsn1[k]
  state.Ts[5] = var.Tsn2[k]

  state.Qhc_o = var.Qhc_o[k]
  set!(state.m_water, m_water)
  set!(state.m_snow, m_snow)
  
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

  mid_res.z_water = z_water
  mid_res.z_snow = z_snow
  mid_res.ρ_snow = ρ_snow[]

  # total GPP . gC/m2/step
  mid_res.GPP = (GPP.o_sunlit + GPP.o_shaded + GPP.u_sunlit + GPP.u_shaded) * 12 * step * 0.000001
  nothing
end
