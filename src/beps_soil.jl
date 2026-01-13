"""
    soil_sm(meteo; kwargs...) -> (results::DataFrame, soil::Soil)

土壤水分与热量独立模拟模块（Standalone Soil Moisture and Heat Simulation）

适用于裸土或稀疏植被（LAI < 2）的土壤过程模拟，不依赖光合作用模块。

# Arguments
- `meteo::DataFrame`: 气象驱动数据，必须包含以下列：
  - `day::Int` - 年积日
  - `hour::Int` - 小时（0-23）
  - `temp::Float64` - 气温 [°C]
  - `rh::Float64` - 相对湿度 [%]
  - `rain::Float64` - 降水 [mm/hour]
  - `Srad::Float64` - 太阳辐射 [W/m²]
  - `wind::Float64` - 风速 [m/s]

# Keyword Arguments
- `SoilType::Int = 8`: 土壤质地类型 (1-11, 见 Init_Soil_Parameters)
- `VegType::Int = 25`: 植被类型（影响 ψ_min, alpha）
- `Tsoil0::Float64 = 10.0`: 初始土壤温度 [°C]
- `θ0::Float64 = 0.35`: 初始土壤湿度 [m³/m³]
- `z_snow0::Float64 = 0.0`: 初始积雪深度 [m]
- `r_drainage::Float64 = 0.5`: 地表排水系数 [0-1]
- `r_root_decay::Float64 = 0.95`: 根系分布衰减因子
- `Tair_annual_mean::Float64 = 10.0`: 年平均气温 [°C]，用于底部边界条件
- `LAI::Float64 = 0.0`: 叶面积指数，用于蒸腾估算
- `canopy_height::Float64 = 0.1`: 冠层高度 [m]，用于空气动力学计算
- `kstep::Float64 = 360.0`: 子时间步长 [秒]

# Returns
- `results::DataFrame`: 土壤状态变量时间序列
- `soil::Soil`: 最终土壤状态对象

# Physics
- Campbell (1974) 5层土壤水分模型，自适应时间步长
- 5层热传导模型，包含冻融动力学
- 完整地表能量平衡（复用 aerodynamic_conductance_jl，含稳定度修正）
- 土壤蒸发 + 简化潜在蒸腾（基于LAI，受土壤水分胁迫调控）

# Example
```julia
using CSV, DataFrames, BEPS

# 读取气象数据
meteo = CSV.read("forcing.csv", DataFrame)

# 运行模拟（裸土）
results, final_state = soil_sm(meteo; SoilType=8, Tsoil0=15.0, θ0=0.35)

# 运行模拟（稀疏植被）
results, final_state = soil_sm(meteo; SoilType=4, LAI=1.0, canopy_height=0.5)
```

# References
- Campbell (1974) - Soil water retention curve
- Chen (2007) Ecological Modelling - Volumetric heat capacity (Eq. 18)
- Bonan (2019) Table 8.2 - Campbell hydraulic parameters
- Brutsaert (1975) - Longwave radiation parameterization
"""
function soil_sm(
    meteo::DataFrame;
    SoilType::Int = 8,
    VegType::Int = 25,
    Tsoil0::Float64 = 10.0,
    θ0::Float64 = 0.35,
    z_snow0::Float64 = 0.0,
    r_drainage::Float64 = 0.5,
    r_root_decay::Float64 = 0.95,
    Tair_annual_mean::Float64 = 10.0,
    LAI::Float64 = 0.0,
    canopy_height::Float64 = 0.1,
    kstep::Float64 = 360.0
)
    # ========== 1. 初始化 ==========
    n = size(meteo, 1)
    soil = Soil()

    # 设置土壤参数和初始状态
    Init_Soil_Parameters(VegType, SoilType, r_root_decay, soil)
    Init_Soil_Status(soil, Tsoil0, meteo.temp[1], θ0, z_snow0)
    soil.r_drainage = r_drainage

    # 状态变量
    z_snow = z_snow0
    z_water = 0.0
    ρ_snow = 250.0  # 初始雪密度 [kg/m³]

    # 地表温度（多层积雪情况）
    Ts0 = meteo.temp[1]      # 地表温度
    Tsm0 = meteo.temp[1]     # 土壤表面温度
    Tsn0 = meteo.temp[1]     # 雪表面温度
    Tsn1 = meteo.temp[1]     # 雪层1温度
    Tsn2 = meteo.temp[1]     # 雪层2温度

    # 分配输出数组
    results = allocate_results(n)

    # ========== 2. 时间步进循环 ==========
    for i in 1:n
        # 2.1 提取气象强迫
        Ta = meteo.temp[i]
        RH = meteo.rh[i]
        prcp = meteo.rain[i]
        Rs = meteo.Srad[i]
        wind = meteo.wind[i]

        # 2.2 降水分配（雨/雪）
        r_rain_g, Δz_snow = partition_precipitation(Ta, prcp, z_snow, ρ_snow)
        z_snow += Δz_snow
        soil.r_rain_g = r_rain_g

        # 2.3 地表能量平衡
        # 计算反照率
        α = calc_soil_albedo(soil.θ[1], soil.θ_sat[1], z_snow)

        # 计算净辐射
        Rn_g = calc_net_radiation(Rs, Ta, RH, α)

        # 计算空气动力学导度（含稳定度修正）
        GH_prev = i > 1 ? results.H_sensible[i-1] : 0.0
        Gheat_g, ra_g = calc_aerodynamic_conductance_bare(
            wind, Ta, canopy_height, LAI, z_snow, GH_prev
        )

        # 2.4 蒸发（复用现有函数）
        perc_snow = Layer3(0.0, 0.0, 0.0)
        mass_snow = Layer3(ρ_snow * z_snow, 0.0, 0.0)

        Evap_soil, Evap_water, Evap_snow, z_water, z_snow =
            evaporation_soil_jl(
                Ta, Ts0, RH, Rn_g, Gheat_g, perc_snow,
                z_water, z_snow, ρ_w * z_water, mass_snow,
                ρ_snow, soil.θ[1], soil.θ_sat[1]
            )

        # 2.5 更新土壤热力学性质
        UpdateSoilThermalConductivity(soil)
        Update_Cs(soil)

        # 2.6 地表温度求解（复用现有函数）
        G, Ts0_new, Tm1, Tsm0_new, Tsn0_new, Tsn1_new, Tsn2_new =
            surface_temperature_jl(
                Ta, RH, z_snow, z_water,
                soil.Cs[1], soil.Cs[1], Gheat_g, soil.dz[1],
                ρ_snow, Ta,  # tempL_u = Ta（裸土无林下层）
                Rn_g, Evap_soil, Evap_water, Evap_snow,
                soil.κ[1], perc_snow.g, soil.G[1],
                Ts0, soil.Tsoil_p[2], soil.Tsoil_p[1], Tsm0,
                Tsn0, Tsn1, Tsn2
            )

        # 更新土壤第一层温度和热通量
        soil.Tsoil_c[1] = Tm1
        soil.G[1] = G
        Ts0, Tsm0, Tsn0, Tsn1, Tsn2 = Ts0_new, Tsm0_new, Tsn0_new, Tsn1_new, Tsn2_new

        # 2.7 土壤热通量（更新第2-5层）
        UpdateHeatFlux(soil, Tair_annual_mean, kstep)

        # 2.8 根系水分吸收（土壤蒸发 + 简化蒸腾）
        # 估算潜在蒸腾
        PET = estimate_potential_transpiration(Rn_g, Ta, RH, Gheat_g, LAI)

        # 根据土壤水分胁迫调整实际蒸腾
        soil_water_factor_v2(soil)
        Trans_actual = PET * soil.f_soilwater

        # 按根系分布分配蒸腾到各层
        Root_Water_Uptake(soil, Trans_actual, 0.0, Evap_soil)

        # 2.9 土壤水分更新
        soil.z_water = z_water
        UpdateSoilMoisture(soil, kstep)
        z_water = soil.z_water

        # 2.10 保存结果
        # 土壤湿度
        results.θ_1[i] = soil.θ[1]
        results.θ_2[i] = soil.θ[2]
        results.θ_3[i] = soil.θ[3]
        results.θ_4[i] = soil.θ[4]
        results.θ_5[i] = soil.θ[5]

        # 土壤温度
        results.Tsoil_1[i] = soil.Tsoil_c[1]
        results.Tsoil_2[i] = soil.Tsoil_c[2]
        results.Tsoil_3[i] = soil.Tsoil_c[3]
        results.Tsoil_4[i] = soil.Tsoil_c[4]
        results.Tsoil_5[i] = soil.Tsoil_c[5]

        # 冰比例
        results.ice_ratio_1[i] = soil.ice_ratio[1]
        results.ice_ratio_2[i] = soil.ice_ratio[2]
        results.ice_ratio_3[i] = soil.ice_ratio[3]
        results.ice_ratio_4[i] = soil.ice_ratio[4]
        results.ice_ratio_5[i] = soil.ice_ratio[5]

        # 地表通量
        results.Evap_soil[i] = Evap_soil
        results.G_surface[i] = G
        results.Ts_ground[i] = Ts0
        results.z_water[i] = z_water
        results.z_snow[i] = z_snow

        # 水量平衡分量
        results.infiltration[i] = soil.r_waterflow[1]
        results.runoff[i] = max(0.0, r_rain_g - soil.r_waterflow[1])
        results.Trans_total[i] = Trans_actual
        results.PET[i] = PET
        results.f_soilwater[i] = soil.f_soilwater

        # 能量平衡分量
        λ = cal_lambda(Ta)
        LE_total = (Evap_soil + Evap_water + Evap_snow + Trans_actual) * λ
        H_sensible = Rn_g - LE_total - G
        results.LE_total[i] = LE_total
        results.H_sensible[i] = H_sensible
    end

    return results, soil
end


# ========== 辅助函数 ==========

"""
    calc_soil_albedo(θ, θ_sat, z_snow) -> α

计算地表反照率（Albedo）

# Arguments
- `θ`: 表层土壤湿度 [m³/m³]
- `θ_sat`: 饱和土壤湿度 [m³/m³]
- `z_snow`: 积雪深度 [m]

# Returns
- `α`: 反照率 [0-1]

# Notes
- 有雪时：α = 0.6
- 无雪时：干土 ~0.25，湿土 ~0.15，线性插值
"""
function calc_soil_albedo(θ::Float64, θ_sat::Float64, z_snow::Float64)
    if z_snow > 0.01
        return 0.6  # 雪的反照率
    else
        # 土壤反照率随湿度变化
        α_dry, α_wet = 0.25, 0.15
        return α_wet + (α_dry - α_wet) * (1.0 - θ / θ_sat)
    end
end


"""
    calc_net_radiation(Rs, Ta, RH, α) -> Rn

计算净辐射（Net Radiation）

# Arguments
- `Rs`: 太阳短波辐射 [W/m²]
- `Ta`: 气温 [°C]
- `RH`: 相对湿度 [%]
- `α`: 反照率 [0-1]

# Returns
- `Rn`: 净辐射 [W/m²]

# Notes
使用 Brutsaert (1975) 大气长波辐射参数化方案
"""
function calc_net_radiation(Rs::Float64, Ta::Float64, RH::Float64, α::Float64)
    # 净短波辐射
    Rns = Rs * (1.0 - α)

    # 长波辐射平衡
    ε = 0.95  # 地表发射率
    σ = 5.67e-8  # Stefan-Boltzmann常数 [W m⁻² K⁻⁴]
    T_K = Ta + 273.15

    # 大气发射率（Brutsaert 1975）
    ea = calc_vapor_pressure(Ta, RH)  # [Pa]
    ε_atm = 1.24 * (ea / T_K)^(1.0/7.0)

    # 净长波辐射（向上为负，向下为正）
    Rnl = ε * σ * T_K^4 * (ε_atm - 1.0)

    return Rns + Rnl
end


"""
    calc_vapor_pressure(Ta, RH) -> ea

计算实际水汽压（Actual Vapor Pressure）

# Arguments
- `Ta`: 气温 [°C]
- `RH`: 相对湿度 [%]

# Returns
- `ea`: 实际水汽压 [Pa]

# Notes
使用 Tetens 公式（WMO 2008）
"""
function calc_vapor_pressure(Ta::Float64, RH::Float64)
    # Tetens 公式：饱和水汽压 [kPa]
    es = 0.6108 * exp(17.27 * Ta / (Ta + 237.3))
    ea = es * RH / 100.0
    return ea * 1000.0  # 转换为 [Pa]
end


"""
    partition_precipitation(Ta, prcp, z_snow, ρ_snow) -> (r_rain_g, Δz_snow)

降水分配：雨 vs 雪

# Arguments
- `Ta`: 气温 [°C]
- `prcp`: 降水强度 [mm/hour]
- `z_snow`: 当前积雪深度 [m]
- `ρ_snow`: 雪密度 [kg/m³]

# Returns
- `r_rain_g`: 到达地表的降雨率 [m/s]
- `Δz_snow`: 新增积雪深度 [m]

# Notes
简化方案：Ta < 0°C 为雪，否则为雨
"""
function partition_precipitation(Ta::Float64, prcp::Float64, z_snow::Float64, ρ_snow::Float64)
    if Ta < 0.0
        # 降雪
        Δz_snow = prcp / 1000.0 / ρ_snow * 3600.0  # mm/h -> m
        r_rain_g = 0.0
    else
        # 降雨
        Δz_snow = 0.0
        r_rain_g = prcp / 1000.0 / 3600.0  # mm/h -> m/s
    end

    return r_rain_g, Δz_snow
end


"""
    estimate_potential_transpiration(Rn, Ta, RH, Gheat, LAI) -> PET

估算潜在蒸腾（Potential Transpiration）

# Arguments
- `Rn`: 净辐射 [W/m²]
- `Ta`: 气温 [°C]
- `RH`: 相对湿度 [%]
- `Gheat`: 空气动力学导度 [m/s]
- `LAI`: 叶面积指数 [-]

# Returns
- `PET`: 潜在蒸腾 [kg/m²/s]

# Notes
- 使用 Penman-Monteith 方程
- 简化冠层气孔导度：gs = 0.01 * LAI
- LAI < 0.01 时返回 0（无植被）
"""
function estimate_potential_transpiration(
    Rn::Float64,
    Ta::Float64,
    RH::Float64,
    Gheat::Float64,
    LAI::Float64
)
    if LAI < 0.01
        return 0.0  # 无植被
    end

    # 气象打包
    met = meteo_pack_jl(Ta, RH)
    (; VPD, Δ, γ) = met
    λ = cal_lambda(Ta)
    cp = cal_cp(Ta, RH)  # 比热容 [J kg⁻¹ K⁻¹]

    # 简化的冠层气孔导度（基于 LAI）
    # 典型草地：LAI=1, gs~0.01 m/s
    gs_canopy = 0.01 * LAI  # [m/s]
    ga = Gheat

    # Penman-Monteith 方程
    PET = 1.0 / λ * (Δ * Rn + ρₐ * cp * VPD * ga) /
          (Δ + γ * (1.0 + ga / gs_canopy))

    return max(0.0, PET)  # [kg/m²/s]
end


"""
    calc_aerodynamic_conductance_bare(wind, Ta, canopy_height, LAI, z_snow, GH_prev)
        -> (Gheat_g, ra_g)

计算裸土/稀疏植被的空气动力学导度

# Arguments
- `wind`: 风速 [m/s]
- `Ta`: 气温 [°C]
- `canopy_height`: 冠层高度 [m]
- `LAI`: 叶面积指数 [-]
- `z_snow`: 积雪深度 [m]
- `GH_prev`: 上一时刻感热通量 [W/m²]，用于稳定度修正

# Returns
- `Gheat_g`: 地表热传导度 [m/s]
- `ra_g`: 地表空气动力学阻抗 [s/m]

# Notes
复用 aerodynamic_conductance_jl 函数，包含完整的稳定度修正
"""
function calc_aerodynamic_conductance_bare(
    wind::Float64,
    Ta::Float64,
    canopy_height::Float64,
    LAI::Float64,
    z_snow::Float64,
    GH_prev::Float64
)
    # 对于裸土/稀疏植被的参数设置
    z_wind = 2.0        # 风速观测高度 [m]
    clumping = 1.0      # 无冠层聚集
    lai_o = LAI
    lai_u = 0.0         # 无林下层

    # 调用现有函数（包含完整的稳定度修正）
    ra_o, ra_u, ra_g, Ga_o, Gb_o, Ga_u, Gb_u =
        aerodynamic_conductance_jl(
            canopy_height, 0.01, z_wind, clumping,
            Ta, wind, GH_prev, lai_o, lai_u
        )

    Gheat_g = 1.0 / ra_g
    return Gheat_g, ra_g
end


"""
    allocate_results(n) -> DataFrame

分配输出结果数组

# Arguments
- `n`: 时间步数

# Returns
- `results`: 包含所有输出变量的 DataFrame
"""
function allocate_results(n::Int)
    DataFrame(
        # 土壤湿度 [m³/m³]
        θ_1 = zeros(n), θ_2 = zeros(n), θ_3 = zeros(n),
        θ_4 = zeros(n), θ_5 = zeros(n),

        # 土壤温度 [°C]
        Tsoil_1 = zeros(n), Tsoil_2 = zeros(n), Tsoil_3 = zeros(n),
        Tsoil_4 = zeros(n), Tsoil_5 = zeros(n),

        # 冰比例 [0-1]
        ice_ratio_1 = zeros(n), ice_ratio_2 = zeros(n), ice_ratio_3 = zeros(n),
        ice_ratio_4 = zeros(n), ice_ratio_5 = zeros(n),

        # 地表通量
        Evap_soil = zeros(n),      # 土壤蒸发 [kg/m²/s]
        G_surface = zeros(n),      # 地表热通量 [W/m²]
        Ts_ground = zeros(n),      # 地表温度 [°C]

        # 水量平衡
        infiltration = zeros(n),   # 入渗 [m/s]
        runoff = zeros(n),         # 径流 [m/s]
        z_water = zeros(n),        # 地表积水 [m]
        z_snow = zeros(n),         # 积雪深度 [m]

        # 蒸散发分量
        Trans_total = zeros(n),    # 总蒸腾 [kg/m²/s]
        PET = zeros(n),            # 潜在蒸腾 [kg/m²/s]
        f_soilwater = zeros(n),    # 土壤水分胁迫因子 [0-1]

        # 能量平衡
        LE_total = zeros(n),       # 潜热通量 [W/m²]
        H_sensible = zeros(n)      # 感热通量 [W/m²]
    )
end


# ========== 导出 ==========
export soil_sm
