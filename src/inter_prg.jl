"""
The inter-module function between main program and modules

# Arguments

- `clumping`  : clumping index
- `param`     : parameter array according to land cover types
- `CosZs`     : cosine of solar zenith angle
- `soilp`     : soil coefficients according to land cover types and soil textures
- `mid_res`   : results struct
"""
function inter_prg_jl(jday::Int, hour::Int, CosZs::T, Ra::Radiation, lai::T, Ω::T,
  forcing::Met, ps::ParamBEPS{T}, state::StateBEPS,
  mid_res::Results, mid_ET::OutputET, cache::LeafCache;
  fix_snowpack::Bool=true, kw...) where {T}

  @unpack Cc_new, Cs_old, Cs_new, Ci_old,
  Tc_old, Tc_new, Gs_old, Gc, Gh, Gw, Gww,
  Gs_new, Ac, Ci_new, Rn, Rns, Rnl,
  leleaf, GPP, LAI, PAI = cache

  # ===== 1. 参数提取和计算 =====
  (; α_canopy_vis, α_canopy_nir,
    α_soil_sat, α_soil_dry, z_canopy_o, z_canopy_u, z_wind,
    g0_w, g1_w, VCmax25, N_leaf, slope_Vc) = ps.veg
  θ_vwp = ps.hydraulic.θ_vwp
  θ_sat = ps.hydraulic.θ_sat

  Vcmax_sunlit, Vcmax_shaded = VCmax(lai, Ω, CosZs, VCmax25, N_leaf, slope_Vc)
  lai_o, lai_u, stem_o, stem_u = lai2!(ps.veg, Ω, CosZs, lai, LAI, PAI)

  # ===== 2. 气象变量初始化 =====
  (; Srad, LR, Tair, RH, wind) = forcing
  precip = forcing.rain / step
  met = meteo_pack_jl(Tair, RH) # 变量类型转换

  # ===== 3. 表面状态初始化 =====
  init_leaf_dbl(Tc_old, Tair - 0.5)

  # 雪水状态 [kg/m² or m]
  m_snow_pre, m_water_pre = state.m_snow, state.m_water
  m_snow, m_water = Layer3(0.0), Layer2()
  z_snow = state.z_snow
  z_water = state.z_water < 0.001 ? 0.0 : state.z_water

  # 雪覆盖和反照率 [-]
  frac_snow, A_snow, frac_water = Layer3(0.0), Layer2(), Layer2()
  α_v = Srad <= 0 ? Layer3() : Layer3(α_canopy_vis)
  α_n = Srad <= 0 ? Layer3() : Layer3(α_canopy_nir)
  ρ_snow, α_v_sw, α_n_sw = init_dbl(state.ρ_snow), init_dbl(), init_dbl()
  Tc = Layer3()

  # 土壤临时变量和中间变量
  radiation_o = radiation_u = radiation_g = ra_g = 0.0

  # ET 相关局部变量
  Trans_o, Trans_u = 0.0, 0.0
  Eil_o, Eil_u = 0.0, 0.0
  EiS_o, EiS_u = 0.0, 0.0
  Evap_soil, Evap_SW, Evap_SS = 0.0, 0.0, 0.0
  Qhc_o, Qhc_u, Qhg = state.Qhc_o, 0.0, 0.0  # Qhc_o 从 state 初始化
  r_rain_g = 0.0

  pai_o = lai_o + stem_o
  pai_u = lai_u + stem_u

  geo_params = (; z_canopy_o, z_canopy_u, z_wind, Ω, lai_o, lai_u, pai_o, pai_u)
  biophys_params = (; g0_w, g1_w, Vcmax_sunlit, Vcmax_shaded)

  # 对雪面温度进行限制
  prev = state.Tsnow_p # 基于地址的修改
  curr = state.Tsnow_c
  clamp!(prev, curr, Tair)
  clamp!(curr, curr, Tair)

  # ===== 4. 亚小时循环 (10步/小时, 360秒/步) =====
  @inbounds for k = 2:kloop+1
    k > 2 && (prev .= curr) # 更新 prev 为上一子时间步的值（k≥3时）

    !fix_snowpack && (ρ_snow[] = 0.0) # TODO: exact as C
    α_v_sw[], α_n_sw[] = 0.0, 0.0

    # /*****  Snowpack stage 1 by X. Luo  *****/
    z_snow = snowpack_stage1_jl(Tair, precip, lai_o, lai_u, Ω,
      m_snow_pre, m_snow, frac_snow, A_snow,
      z_snow, ρ_snow, α_v_sw, α_n_sw)

    # /*****  Rainfall stage 1 by X. Luo  *****/
    r_rain_g = rainfall_stage1_jl(Tair, precip, frac_water, m_water, m_water_pre, lai_o, lai_u, Ω)

    # 土壤反照率计算 [-]
    α_g = if state.θ_prev[2] < θ_vwp[2] * 0.5
      α_soil_dry
    else
      (state.θ_prev[2] - θ_vwp[2] * 0.5) / (θ_sat[2] - θ_vwp[2] * 0.5) *
      (α_soil_sat - α_soil_dry) + α_soil_dry
    end
    α_v.g = 2.0 / 3.0 * α_g
    α_n.g = 4.0 / 3.0 * α_g

    # /*****  Soil water factor module by L. He  *****/
    soil_water_factor_v2(state, ps)
    f_soilwater = min(state.f_soilwater, 1.0) # used in `photosynthesis`

    # 感热通量初值用于空气动力学导度计算 [W/m²]
    H_canopy_o = Qhc_o  # 使用上一步的值

    init_leaf_dbl(Ci_old, 0.7 * CO2_air)
    init_leaf_dbl2(Gs_old, 1.0 / 200.0, 1.0 / 300.0)

    perc_snow_o = A_snow.o / lai_o / 2 # 上层冠层雪覆盖分数
    perc_snow_u = A_snow.u / lai_u / 2 # 下层冠层雪覆盖分数

    Tc.g = Tair   # 地表温度初值用气温代替

    # 能量平衡迭代求解冠层温度
    # /*****  Canopy Energy Balance Iteration extracted by AI Agent  *****/
    snow_params = (; perc_snow_o, perc_snow_u, frac_snow, frac_water, α_v_sw, α_n_sw, α_v, α_n)

    radiation_o, radiation_u, radiation_g, ra_g, _ = solve_canopy_energy_balance!(
      cache, met, forcing, geo_params, snow_params, biophys_params,
      Tc, Ra, H_canopy_o, CosZs, f_soilwater, frac_water
    )
    multiply!(GPP, Ac, LAI)

    Trans_o, Trans_u = transpiration_jl(Tc_new, Tair, RH, Gw, LAI) # X. Luo

    # /*****  Evaporation and sublimation from canopy by X. Luo  *****/
    Eil_o, Eil_u, EiS_o, EiS_u = evaporation_canopy_jl(Tc_new, Tair, RH,
      Gww, PAI, frac_water, frac_snow)

    rainfall_stage2_jl(Eil_o, Eil_u, m_water) # X. Luo
    set!(m_water_pre, m_water)

    snowpack_stage2_jl(EiS_o, EiS_u, m_snow) # X. Luo

    # /*****  Evaporation from soil module by X. Luo  *****/
    # ra_g 是地表到参考高度的总阻抗 (地表→下层冠层→上层冠层→参考高度)
    Gheat_g = 1 / ra_g  # 地表空气动力学传热导度 [m/s]
    mass_water_g = ρ_w * z_water  # 地表水质量 [kg/m²]

    Evap_soil, Evap_SW, Evap_SS, z_water, z_snow =
      evaporation_soil_jl(Tair, prev.T_surf, RH, radiation_g, Gheat_g,
        frac_snow, z_water, z_snow, mass_water_g, m_snow,
        ρ_snow[], state.θ_prev[1], θ_sat[1])

    # /*****  Surface temperature by X. Luo  *****/
    state.G[1] = surface_temperature!(state, ps, prev, curr,
      radiation_g, Tc.u, Tair, RH, z_snow, z_water,
      ρ_snow[], frac_snow.g, Gheat_g,
      Evap_soil, Evap_SW, Evap_SS)

    # /*****  Snowpack stage 3 by X. Luo  *****/
    z_snow, z_water = snowpack_stage3_jl(Tair, curr.T_snow0, prev.T_snow0,
      ρ_snow[], z_snow, z_water, m_snow)
    set!(m_snow_pre, m_snow)
    state.z_snow = z_snow

    # /*****  Sensible heat flux by X. Luo  *****/
    Qhc_o, Qhc_u, Qhg = sensible_heat_jl(Tc_new, curr.T_surf, Tair, RH, Gh, Gheat_g, PAI)
    UpdateHeatFlux(state, Tair, kstep)

    # /*****  Soil water module by L. He  *****/
    Root_Water_Uptake(state, Trans_o, Trans_u, Evap_soil)

    state.r_rain_g = r_rain_g
    state.z_water = z_water

    UpdateSoilMoisture(state, ps, kstep)
    z_water = state.z_water
  end  # end of sub-hourly loop

  # ===== 5. 时间步结束：状态更新 =====
  state.Qhc_o = Qhc_o
  set!(state.m_water, m_water)
  set!(state.m_snow, m_snow)
  state.ρ_snow = ρ_snow[]

  # ===== 6. 输出结果汇总 =====
  @pack! mid_ET = Trans_o, Trans_u, Eil_o, Eil_u, EiS_o, EiS_u,
  Evap_soil, Evap_SW, Evap_SS, Qhc_o, Qhc_u, Qhg
  update_ET!(mid_ET, mid_res, Tair)

  mid_res.Net_Rad = radiation_o + radiation_u + radiation_g
  mid_res.gpp_o_sunlit = GPP.o_sunlit
  mid_res.gpp_u_sunlit = GPP.u_sunlit
  mid_res.gpp_o_shaded = GPP.o_shaded
  mid_res.gpp_u_shaded = GPP.u_shaded

  mid_res.z_water = z_water
  mid_res.z_snow = z_snow
  mid_res.ρ_snow = ρ_snow[]

  GPP = GPP.o_sunlit + GPP.o_shaded + GPP.u_sunlit + GPP.u_shaded
  mid_res.GPP = GPP * 12 * step * 1e-6  # [umol m-2 s-1] -> [gC m-2]
  nothing
end


"""
    solve_canopy_energy_balance!(cache, met, forcing, geo_params, snow_params, biophys_params,
        Tc, Ra, H_canopy_o, CosZs, f_water)

Iteratively solves the canopy energy balance to determine canopy temperatures and fluxes.
Extracted from `inter_prg_jl` to improve readability and maintainability.
"""
function solve_canopy_energy_balance!(
  cache::LeafCache, met::NamedTuple, forcing::Met,
  geo_params, snow_params, biophys_params,
  Tc, Ra::Radiation, H_canopy_o::Float64, CosZs::T, f_soilwater::T, frac_water::Layer2{T}
) where {T}

  # Unpack required variables
  @unpack Cc_new, Cs_old, Cs_new, Ci_old,
  Tc_old, Tc_new, Gs_old, Gc, Gh, Gw, Gww,
  Gs_new, Ac, Ci_new, Rn, Rns, Rnl,
  leleaf, PAI = cache

  (; ρₐ, cp, VPD, ea, Δ, γ) = met
  (; Tair, RH, wind, LR) = forcing
  (; z_canopy_o, z_canopy_u, z_wind, Ω, lai_o, lai_u, pai_o, pai_u) = geo_params
  (; perc_snow_o, perc_snow_u, frac_snow, α_v_sw, α_n_sw, α_v, α_n) = snow_params
  (; g0_w, g1_w, Vcmax_sunlit, Vcmax_shaded) = biophys_params

  Srad = forcing.Srad # Directly use Srad

  radiation_o = radiation_u = radiation_g = ra_g = 0.0
  n_iter = 0

  while true
    n_iter += 1
    # /***** Aerodynamic conductance module by G.Mo  *****/
    ra_o, ra_u, ra_g, Ga_o, Gb_o, Ga_u, Gb_u =
      aerodynamic_conductance_jl(z_canopy_o, z_canopy_u, z_wind, Ω, Tair, wind, H_canopy_o,
        pai_o, pai_u)

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
      lai_o, lai_u, pai_o, pai_u, PAI,
      Ω, Tair, RH, LR,
      α_v_sw[], α_n_sw[], α_v, α_n,
      perc_snow_o, perc_snow_u, frac_snow.g,
      Rn, Rns, Rnl, Ra)

    # /*****  Photosynthesis module by B. Chen  *****/
    update_Gw!(Gw, Gs_old, Ga_o, Ga_u, Gb_o, Gb_u) # 水汽导度
    latent_heat!(leleaf, Gw, VPD, Δ, Tc_old, Tair, ρₐ, cp, γ)

    if (CosZs > 0)
      photosynthesis(Tc_old, Rns, Ci_old, leleaf,
        Tair, ea, f_soilwater, g0_w, g1_w,
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
    Leaf_Temperatures_jl(Tair, Δ, γ, VPD, cp,
      Gw, Gww, Gh, frac_water, frac_snow, Rn, Tc_new)

    # 计算上层冠层感热通量用于下次迭代 [W/m²]
    H_o_sunlit = (Tc_new.o_sunlit - Tair) * ρₐ * cp * Gh.o_sunlit
    H_o_shaded = (Tc_new.o_shaded - Tair) * ρₐ * cp * Gh.o_shaded
    H_canopy_o = H_o_sunlit * PAI.o_sunlit + H_o_shaded * PAI.o_shaded

    # 检查冠层温度是否收敛 (精度0.02°C)
    if (abs(Tc_new.o_sunlit - Tc_old.o_sunlit) < 0.02 &&
        abs(Tc_new.o_shaded - Tc_old.o_shaded) < 0.02 &&
        abs(Tc_new.u_sunlit - Tc_old.u_sunlit) < 0.02 &&
        abs(Tc_new.u_shaded - Tc_old.u_shaded) < 0.02)
      break
    else
      if (n_iter > 22)  # 迭代未收敛，使用气温作为冠层温度
        init_leaf_dbl(Tc_old, Tair)
        break
      else
        set!(Tc_old, Tc_new)
      end
    end
  end

  return radiation_o, radiation_u, radiation_g, ra_g, H_canopy_o
end
