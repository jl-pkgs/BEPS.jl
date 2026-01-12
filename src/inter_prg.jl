"""
The inter-module function between main program and modules

# Arguments

- `jday`      : day of year
- `hour`      : hour of day
- `lai`       : leaf area index
- `clumping`  : clumping index
- `param`     : parameter array according to land cover types
- `meteo`     : meteorological data
- `CosZs`     : cosine of solar zenith angle
- `var_o`     : temporary variables array of last time step
- `var_n`     : temporary variables array of this time step
- `soilp`     : soil coefficients according to land cover types and soil textures
- `mid_res`   : results struct
"""
function inter_prg_jl(
  jday::Int, hour::Int,
  lai::T, Ω::T, param::VegParam{T}, meteo::Met, CosZs::T,
  state::State{T}, soil::AbstractSoil,
  Ra::Radiation,
  mid_res::Results, mid_ET::OutputET, var::InterTempVars; 
  VegType::Int=1, 
  fix_snowpack::Bool=true, kw...) where {T}

  init_vars!(var)
  @unpack Cc_new, Cs_old, Cs_new, Ci_old,
  Tc_old, Tc_new, Gs_old, Gc, Gh, Gw, Gww,
  Gs_new, Ac, Ci_new, Rn, Rns, Rnl,
  leleaf, GPP, LAI, PAI = var.TempLeafs

  # 土壤分层临时变量
  dz = zeros(layer + 1)
  κ = zeros(layer + 2)

  # 空气动力学阻抗和导度 [s/m] 或 [mol/m²/s]
  ra_g = 0.0
  Ga_o = 0.0
  Gb_o = 0.0
  Ga_u = 0.0
  Gb_u = 0.0

  # 净辐射 [W/m²]
  radiation_o, radiation_u, radiation_g = 0.0, 0.0, 0.0

  # /*****  Vcmax-Nitrogen calculations by G.Mo  *****/
  (; LAI_max_o, LAI_max_u,
    α_canopy_vis, α_canopy_nir,
    α_soil_sat, α_soil_dry,
    z_canopy_o, z_canopy_u, z_wind,
    g0_w, g1_w,
    VCmax25, N_leaf, slope_Vc) = param

  Vcmax_sunlit, Vcmax_shaded = VCmax(lai, Ω, CosZs, VCmax25, N_leaf, slope_Vc)

  # /*****  LAI calculation module by B. Chen  *****/
  lai_o = lai < 0.1 ? 0.1 : lai
  if (VegType == 25 || VegType == 40) # C4植物参数
    lai_u = 0.01
  else
    lai_u = 1.18 * exp(-0.99 * lai_o)
  end
  (lai_u > lai_o) && (lai_u = 0.01)

  stem_o = LAI_max_o * 0.2    # 上层茎面积指数
  stem_u = LAI_max_u * 0.2    # 下层茎面积指数

  # 分离阳叶和阴叶LAI
  lai2(Ω, CosZs, stem_o, stem_u, lai_o, lai_u, LAI, PAI)

  # /*****  Initialization of this time step  *****/
  Srad = meteo.Srad      # 入射短波辐射 [W/m²]
  RH = meteo.rh          # 相对湿度 [0-1]
  wind = meteo.wind      # 风速 [m/s]
  precip = meteo.rain / step  # 降水速率 [m/s]
  T_air = meteo.temp     # 气温 [°C]

  met = meteo_pack_jl(T_air, RH)
  (; Δ, γ, cp, VPD, ea) = met

  init_leaf_dbl(Tc_old, T_air - 0.5)

  # Ground surface temperature
  # [Ts0, Tsn, Tsm0, Tsn1, Tsn2]
  var.Ts0[1] = clamp(state.Ts[1], T_air - 2.0, T_air + 2.0)  # ground0
  var.Tsn0[1] = clamp(state.Ts[2], T_air - 2.0, T_air + 2.0) # snow0
  var.Tsm0[1] = clamp(state.Ts[3], T_air - 2.0, T_air + 2.0) # any
  var.Tsn1[1] = clamp(state.Ts[4], T_air - 2.0, T_air + 2.0) # snow1
  var.Tsn2[1] = clamp(state.Ts[5], T_air - 2.0, T_air + 2.0) # snow2
  var.Qhc_o[1] = state.Qhc_o

  # 雪和水的质量状态 [kg/m²]
  m_snow_pre = state.m_snow
  m_water_pre = state.m_water
  m_snow = Layer3(0.0)
  m_water = Layer2()

  # 雪和水的深度 [m]
  z_snow = soil.z_snow
  z_water = soil.z_water
  (z_water < 0.001) && (z_water = 0.0)

  # 雪覆盖分数 [-]
  f_snow = Layer3(0.0)
  A_snow = Layer2()
  f_water = Layer2()

  # 反照率 [-]
  if (Srad <= 0)
    α_v = Layer3()
    α_n = Layer3()
  else
    α_v = Layer3(α_canopy_vis)
    α_n = Layer3(α_canopy_nir)
  end

  # 雪密度和反照率
  ρ_snow = init_dbl(state.ρ_snow)
  α_v_sw = init_dbl()
  α_n_sw = init_dbl()

  # 冠层温度 [°C]
  Tc = Layer3()

  # 十次亚小时循环: 每小时分10步, 每步360秒
  @inbounds for k_step = 2:kloop+1
    !fix_snowpack && (ρ_snow[] = 0.0) # TODO: exact as C
    α_v_sw[], α_n_sw[] = 0.0, 0.0

    # /*****  Snowpack stage 1 by X. Luo  *****/
    z_snow = snowpack_stage1_jl(T_air, precip,
      lai_o, lai_u, Ω,
      m_snow_pre, m_snow, f_snow, A_snow,
      z_snow, ρ_snow,
      α_v_sw, α_n_sw)

    # /*****  Rainfall stage 1 by X. Luo  *****/
    var.r_rain_g[k_step] = rainfall_stage1_jl(T_air, precip, f_water, m_water, m_water_pre, lai_o, lai_u, Ω)

    if (soil.θ_prev[2] < soil.θ_vwp[2] * 0.5)
      α_g = α_soil_dry
    else
      α_g = (soil.θ_prev[2] - soil.θ_vwp[2] * 0.5) / (soil.θ_sat[2] - soil.θ_vwp[2] * 0.5) * (α_soil_sat - α_soil_dry) + α_soil_dry
    end

    α_v.g = 2.0 / 3.0 * α_g
    α_n.g = 4.0 / 3.0 * α_g

    # /*****  Soil water factor module by L. He  *****/
    soil_water_factor_v2(soil)
    f_soilwater = min(soil.f_soilwater, 1.0) # used in `photosynthesis`

    # 感热通量初值用于空气动力学导度计算 [W/m²]
    H_canopy_o = var.Qhc_o[k_step-1]

    init_leaf_dbl(Ci_old, 0.7 * CO2_air)
    init_leaf_dbl2(Gs_old, 1.0 / 200.0, 1.0 / 300.0)

    perc_snow_o = A_snow.o / lai_o / 2 # 上层冠层雪覆盖分数
    perc_snow_u = A_snow.u / lai_u / 2 # 下层冠层雪覆盖分数

    Tc.g = T_air   # 地表温度初值用气温代替

    # 能量平衡迭代求解冠层温度
    n_iter = 0
    while true
      n_iter += 1
      # /***** Aerodynamic conductance module by G.Mo  *****/
      ra_o, ra_u, ra_g, Ga_o, Gb_o, Ga_u, Gb_u =
        aerodynamic_conductance_jl(z_canopy_o, z_canopy_u, z_wind, Ω, T_air, wind, H_canopy_o,
          lai_o + stem_o, lai_u + stem_u)

      # 热量传输导度 [mol/m²/s]
      init_leaf_dbl2(Gh,
        1.0 / (1.0 / Ga_o + 0.5 / Gb_o),
        1.0 / (1.0 / Ga_u + 0.5 / Gb_u))
      # 水汽传输导度 (湿表面) [mol/m²/s]
      init_leaf_dbl2(Gww,
        1.0 / (1.0 / Ga_o + 1.0 / Gb_o + 100),
        1.0 / (1.0 / Ga_u + 1.0 / Gb_u + 100))

      # 上下层冠层平均温度 [°C]
      Tc.o = (Tc_old.o_sunlit * PAI.o_sunlit + Tc_old.o_shaded * PAI.o_shaded) /
             (PAI.o_sunlit + PAI.o_shaded)
      Tc.u = (Tc_old.u_sunlit * PAI.u_sunlit + Tc_old.u_shaded * PAI.u_shaded) /
             (PAI.u_sunlit + PAI.u_shaded)

      # /*****  Net radiation at canopy and leaf level module by X.Luo  *****/
      radiation_o, radiation_u, radiation_g = netRadiation_jl(
        Srad, CosZs, Tc,
        lai_o, lai_u, lai_o + stem_o, lai_u + stem_u, PAI,
        Ω, T_air, RH,
        α_v_sw[], α_n_sw[], α_v, α_n,
        perc_snow_o, perc_snow_u, f_snow.g,
        Rn, Rns, Rnl, Ra)

      # /*****  Photosynthesis module by B. Chen  *****/
      update_Gw!(Gw, Gs_old, Ga_o, Ga_u, Gb_o, Gb_u) # 水汽导度
      latent_heat!(leleaf, Gw, VPD, Δ, Tc_old, T_air, ρₐ, cp, γ)

      if (CosZs > 0)
        photosynthesis(Tc_old, Rns, Ci_old, leleaf,
          T_air, ea, f_soilwater, g0_w, g1_w,
          Gb_o, Gb_u, Vcmax_sunlit, Vcmax_shaded,
          Gs_new, Ac, Ci_new; version="julia")
      else
        init_leaf_dbl(Gs_new, 0.0001)
        init_leaf_dbl(Ac, 0.0)
        init_leaf_dbl(Ci_new, CO2_air * 0.7)

        init_leaf_dbl(Cs_new, CO2_air)
        init_leaf_dbl(Cc_new, CO2_air * 0.7 * 0.8)
      end

      set!(Ci_old, Ci_new)
      set!(Cs_old, Cs_new)
      set!(Gs_old, Gs_new)

      update_Gw!(Gw, Gs_new, Ga_o, Ga_u, Gb_o, Gb_u)
      update_Gc!(Gc, Gs_new, Ga_o, Ga_u, Gb_o, Gb_u)

      # /***** Leaf temperatures module by L. He  *****/
      Leaf_Temperatures_jl(T_air, Δ, γ, VPD, cp,
        Gw, Gww, Gh, f_water, f_snow, Rn, Tc_new)

      # 计算上层冠层感热通量用于下次迭代 [W/m²]
      H_o_sunlit = (Tc_new.o_sunlit - T_air) * ρₐ * cp * Gh.o_sunlit
      H_o_shaded = (Tc_new.o_shaded - T_air) * ρₐ * cp * Gh.o_shaded
      H_canopy_o = H_o_sunlit * PAI.o_sunlit + H_o_shaded * PAI.o_shaded

      # 检查冠层温度是否收敛 (精度0.02°C)
      if (abs(Tc_new.o_sunlit - Tc_old.o_sunlit) < 0.02 &&
          abs(Tc_new.o_shaded - Tc_old.o_shaded) < 0.02 &&
          abs(Tc_new.u_sunlit - Tc_old.u_sunlit) < 0.02 &&
          abs(Tc_new.u_shaded - Tc_old.u_shaded) < 0.02)
        break
      else
        if (n_iter > 22)  # 迭代未收敛，使用气温作为冠层温度
          init_leaf_dbl(Tc_old, T_air)
          break
        else
          set!(Tc_old, Tc_new)
        end
      end
    end  # end of energy balance iteration

    multiply!(GPP, Ac, LAI)

    # /*****  Transpiration by X. Luo  *****/
    var.Trans_o[k_step], var.Trans_u[k_step] = transpiration_jl(Tc_new, T_air, RH, Gw, LAI)

    # /*****  Evaporation and sublimation from canopy by X. Luo  *****/
    var.Eil_o[k_step], var.Eil_u[k_step], var.EiS_o[k_step], var.EiS_u[k_step] = evaporation_canopy_jl(
      Tc_new, T_air, RH,
      Gww, PAI, f_water, f_snow)

    # /*****  Rainfall stage 2 by X. Luo  *****/
    rainfall_stage2_jl(var.Eil_o[k_step], var.Eil_u[k_step], m_water)
    set!(m_water_pre, m_water)

    # /*****  Snowpack stage 2 by X. Luo  *****/
    snowpack_stage2_jl(var.EiS_o[k_step], var.EiS_u[k_step], m_snow)

    # /*****  Evaporation from soil module by X. Luo  *****/
    Gheat_g = 1 / ra_g  # 地表传热导度 [mol/m²/s]
    mass_water_g = ρ_w * z_water  # 地表水质量 [kg/m²]

    var.Evap_soil[k_step], var.Evap_SW[k_step], var.Evap_SS[k_step], z_water, z_snow =
      evaporation_soil_jl(Tc.g, var.Ts0[k_step-1], RH, radiation_g, Gheat_g,
        f_snow,
        z_water, z_snow, mass_water_g, m_snow,
        ρ_snow[], soil.θ_prev[1], soil.θ_sat[1])

    # /*****  Soil Thermal Conductivity module by L. He  *****/
    UpdateSoilThermalConductivity(soil)
    Update_Cs(soil)

    # /*****  Surface temperature by X. Luo  *****/
    var.Cs .= 0.0
    var.Tm .= 0.0
    var.G .= 0.0

    var.Cs[1, k_step] = soil.Cs[1]
    var.Cs[2, k_step] = soil.Cs[1]
    var.Tc_u[k_step] = Tc.u
    κ[2] = soil.κ[1]
    dz[2] = soil.dz[1]

    var.Tm[1, k_step-1] = soil.Tsoil_p[1]
    var.Tm[2, k_step-1] = soil.Tsoil_p[2]
    var.G[2, k_step] = soil.G[1]

    var.G[1, k_step], var.Ts0[k_step], var.Tm[1, k_step], var.Tsm0[k_step],
    var.Tsn0[k_step], var.Tsn1[k_step], var.Tsn2[k_step] =
      surface_temperature_jl(T_air, RH, z_snow, z_water,
        var.Cs[2, k_step], var.Cs[1, k_step], Gheat_g, dz[2], ρ_snow[], var.Tc_u[k_step],
        radiation_g, var.Evap_soil[k_step], var.Evap_SW[k_step], var.Evap_SS[k_step],
        κ[2],
        f_snow.g, var.G[2, k_step],
        var.Ts0[k_step-1],
        var.Tm[2, k_step-1], var.Tm[1, k_step-1], var.Tsm0[k_step-1],
        var.Tsn0[k_step-1], var.Tsn1[k_step-1], var.Tsn2[k_step-1])

    soil.Tsoil_c[1] = var.Tm[1, k_step]

    # /*****  Snowpack stage 3 by X. Luo  *****/
    z_snow, z_water = snowpack_stage3_jl(T_air, var.Tsn0[k_step], var.Tsn0[k_step-1],
      ρ_snow[], z_snow, z_water, m_snow)
    set!(m_snow_pre, m_snow)

    # /*****  Sensible heat flux by X. Luo  *****/
    var.Qhc_o[k_step], var.Qhc_u[k_step], var.Qhg[k_step] =
      sensible_heat_jl(Tc_new, var.Ts0[k_step], T_air, RH,
        Gh, Gheat_g, PAI)

    # /*****  Soil water module by L. He  *****/
    soil.z_snow = z_snow
    soil.G[1] = var.G[1, k_step]

    UpdateHeatFlux(soil, T_air, kstep)
    Root_Water_Uptake(soil, var.Trans_o[k_step], var.Trans_u[k_step], var.Evap_soil[k_step])

    soil.r_rain_g = var.r_rain_g[k_step]
    soil.z_water = z_water

    UpdateSoilMoisture(soil, kstep)
    z_water = soil.z_water
  end  # end of sub-hourly loop

  k_step = kloop + 1  # 最后一步索引
  var.Tsn1[k_step] = clamp(var.Tsn1[k_step], -40.0, 40.0)
  var.Tsn2[k_step] = clamp(var.Tsn2[k_step], -40.0, 40.0)

  # 更新状态变量
  state.Ts[1] = var.Ts0[k_step]
  state.Ts[2] = var.Tsn0[k_step]
  state.Ts[3] = var.Tsm0[k_step]
  state.Ts[4] = var.Tsn1[k_step]
  state.Ts[5] = var.Tsn2[k_step]

  state.Qhc_o = var.Qhc_o[k_step]
  set!(state.m_water, m_water)
  set!(state.m_snow, m_snow)
  state.ρ_snow = ρ_snow[]

  # 输出结果汇总
  mid_res.Net_Rad = radiation_o + radiation_u + radiation_g

  OutputET!(mid_ET,
    var.Trans_o, var.Trans_u,
    var.Eil_o, var.Eil_u,
    var.EiS_o, var.EiS_u,
    var.Evap_soil, var.Evap_SW, var.Evap_SS, var.Qhc_o, var.Qhc_u, var.Qhg, k_step)
  update_ET!(mid_ET, mid_res, T_air)

  mid_res.gpp_o_sunlit = GPP.o_sunlit   # [umol C/m2/s], !note
  mid_res.gpp_u_sunlit = GPP.u_sunlit
  mid_res.gpp_o_shaded = GPP.o_shaded
  mid_res.gpp_u_shaded = GPP.u_shaded

  mid_res.z_water = z_water
  mid_res.z_snow = z_snow
  mid_res.ρ_snow = ρ_snow[]

  GPP = GPP.o_sunlit + GPP.o_shaded + GPP.u_sunlit + GPP.u_shaded
  mid_res.GPP = GPP * 12 * step * 1e-6 # [umol m-2 s-1] -> [gC m-2]
  nothing
end
