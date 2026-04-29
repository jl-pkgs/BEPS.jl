using BEPS
using Test
using DataFrames

@testset "独立光合模块测试" begin
  @testset "参数初始化" begin
    params = InitParam_Photo_Farquhar()

    @test params.Vcmax25 > 0
    @test params.Jmax25 > 0
    @test params.evc > 0
    @test params.g0_w >= 0
    @test params.g1_w > 0
    @test params.clumping > 0
    # ko25 单位：mmol/mol，应与 O2 浓度（210 mmol/mol）同量级
    @test params.ko25 > 100 && params.ko25 < 1000
    # tau25 为 Rubisco 特异性因子，典型范围 2000-3000
    @test params.tau25 > 1000 && params.tau25 < 5000
  end

  @testset "基本功能" begin
    params = InitParam_Photo_Farquhar()

    result = BEPS.Photosynthesis.photosynthesis(25.0, 60.0, 500.0, 2.0, params)

    @test result.An.sunlit > 0
    @test result.An.shaded >= 0
    @test result.Gs.sunlit > 0
    @test result.Gs.shaded >= 0
    # 高光（Srad=500 W/m²）时阳/阴叶都受 Rubisco 限制，An 相同
    @test result.An.sunlit >= result.An.shaded
    # An 应在物理合理范围内：C3 植物典型值 2–30 μmol m-2 s-1
    @test 1.0 < result.An.sunlit < 40.0
    # 胞间 CO2 应为正值且低于大气 CO2
    @test 0.0 < result.Ci.sunlit < 380.0
    # 气孔导度应在合理范围：0.001–1.0 mol m-2 s-1
    @test 0.0 < result.Gs.sunlit < 1.0
  end

  @testset "弱光下阳生 > 阴生" begin
    params = InitParam_Photo_Farquhar()

    # 弱光（Srad=100 W/m²）时阴叶受光限制，阳叶 An 高于阴叶
    result_low = BEPS.Photosynthesis.photosynthesis(25.0, 60.0, 100.0, 2.0, params)
    @test result_low.An.sunlit > result_low.An.shaded
  end

  @testset "温度响应" begin
    params = InitParam_Photo_Farquhar()

    # 本参数化方案（Bernacchi kc25=404，Medlyn 温度响应）最适温度约 15°C
    # 验证 bell-shaped 响应：最适温度附近 An 高于极端高低温
    result_cold = BEPS.Photosynthesis.photosynthesis(-5.0, 60.0, 500.0, 2.0, params)
    result_opt  = BEPS.Photosynthesis.photosynthesis(15.0, 60.0, 500.0, 2.0, params)
    result_hot  = BEPS.Photosynthesis.photosynthesis(40.0, 60.0, 500.0, 2.0, params)

    @test result_opt.An.sunlit > result_cold.An.sunlit   # 低温（-5°C）抑制
    @test result_opt.An.sunlit > result_hot.An.sunlit    # 高温（40°C）抑制
    @test result_cold.An.sunlit > 0                       # 低温下仍有正 An
  end

  @testset "光响应" begin
    params = InitParam_Photo_Farquhar()

    result_dark  = BEPS.Photosynthesis.photosynthesis(25.0, 60.0, 0.0,   2.0, params)
    result_light = BEPS.Photosynthesis.photosynthesis(25.0, 60.0, 500.0, 2.0, params)

    @test result_light.An.sunlit > result_dark.An.sunlit
    @test result_dark.An.sunlit < 0  # 暗呼吸占主导，净光合为负
  end

  @testset "水分胁迫" begin
    params = InitParam_Photo_Farquhar()

    result_no_stress = BEPS.Photosynthesis.photosynthesis(25.0, 60.0, 500.0, 2.0, params; β_soil=1.0)
    result_stress    = BEPS.Photosynthesis.photosynthesis(25.0, 60.0, 500.0, 2.0, params; β_soil=0.3)

    @test result_no_stress.An.sunlit > result_stress.An.sunlit
  end

  @testset "物理正确性" begin
    params = InitParam_Photo_Farquhar()
    rugc  = 8.314
    TK25  = 298.16
    T25   = 25.0 + 273.15

    # Γstar（CO2 补偿点）：C3 植物 25°C 约 35–45 μmol mol-1
    TBOLTZ_f = (T, E) -> exp((T - TK25) * E / (T * TK25 * rugc))
    tau = params.tau25 * TBOLTZ_f(T25, -29000.0)
    Γstar = 105000.0 / tau
    @test 30.0 < Γstar < 50.0

    # O2/Ko 比例：25°C 约 0.6–1.0
    Ko = params.ko25 * TBOLTZ_f(T25, 14500.0)
    @test 0.5 < (210.0 / Ko) < 1.2

    # Vcmax 和 Jmax 温度响应需在 25°C 归一化
    Hd  = 200000.0; S = 640.0
    fTv_f = (T, Vc, evc) -> begin
      ft  = TBOLTZ_f(T, evc)
      num = 1.0 + exp((S * T25 - Hd) / (T25 * rugc))
      den = 1.0 + exp((S * T  - Hd) / (T  * rugc))
      Vc * ft * num / den
    end
    Vcmax_at_25 = fTv_f(T25, params.Vcmax25, params.evc)
    @test isapprox(Vcmax_at_25, params.Vcmax25, rtol=0.01)  # 1% 以内
  end

  @testset "与原始 photosynthesis_jl 对比" begin
    #=
    对比说明：
    - 原始 photosynthesis_jl：使用 TBOLTZ 钟形温度响应（topt=301K），
      Michaelis 常数 kc25=274.6, ko25=419.8, tau25=2904（Harley & Baldocchi 1995），
      最适温度 ~25–28°C
    - 独立 Photosynthesis 模块：使用 Medlyn et al. 2002 归一化 Arrhenius 温度响应，
      Michaelis 常数 kc25=404, ko25=248, tau25=2600（Bernacchi et al. 2001），
      最适温度 ~15°C（北方针叶林参数化）
    两者 Vcmax25 和 Jmax25 相同，但温度响应函数和 Michaelis 常数不同。
    =#
    params = InitParam_Photo_Farquhar()
    g0_w = params.g0_w; g1_w = params.g1_w
    Vcmax25 = params.Vcmax25

    # gb_w 单位转换：standalone 用 0.5 mol/m²/s，原始用 s/m
    # 0.5 mol/m²/s ≈ 0.5 * T_K * pstat273 s/m，pstat273 = 0.022624/(273.16*1.013)
    pstat273 = 0.022624 / (273.16 * 1.013)
    RH = 60.0; Srad = 500.0

    # 在等价单叶条件下对比（LAI→0 ≈ 单叶，取阳生叶）
    Ts = [0.0, 10.0, 20.0, 25.0, 30.0]
    for T in Ts
      T_K = T + 273.15
      ea  = RH / 100 * 0.6108 * exp(17.27 * T / (T + 237.3))
      gb_w_s = 0.5 * T_K * pstat273           # s/m

      _, An_orig, ci_orig = photosynthesis_jl(
        T_K, Srad, ea, gb_w_s, Vcmax25, 1.0, g0_w, g1_w, 0.7 * 380.0, T, 0.0
      )
      r      = BEPS.Photosynthesis.photosynthesis(T, RH, Srad, 0.01, params)
      An_new = r.An.sunlit
      ci_new = r.Ci.sunlit

      # 两者均应给出正 An（在 0–30°C 范围内）
      @test An_orig >= 0.0
      @test An_new  >= 0.0
      # ci 应在合理范围内
      @test 0.0 < ci_orig < 380.0
      @test 0.0 < ci_new  < 380.0
    end

    # 原始模块在 25°C 最接近最适，应大于在 0°C
    ea_25 = RH / 100 * 0.6108 * exp(17.27 * 25.0 / (25.0 + 237.3))
    ea_0  = RH / 100 * 0.6108 * exp(17.27 * 0.0  / (0.0  + 237.3))
    _, An_orig_25, _ = photosynthesis_jl(298.15, Srad, ea_25,
                                          0.5 * 298.15 * pstat273, Vcmax25, 1.0, g0_w, g1_w,
                                          0.7 * 380.0, 25.0, 0.0)
    _, An_orig_0, _  = photosynthesis_jl(273.15, Srad, ea_0,
                                          0.5 * 273.15 * pstat273, Vcmax25, 1.0, g0_w, g1_w,
                                          0.7 * 380.0, 0.0, 0.0)
    @test An_orig_25 > An_orig_0  # 原始模块：25°C > 0°C（峰值在 ~28°C）
  end

  @testset "参数堆叠" begin
    params = InitParam_Photo_Farquhar()
    df = parameters(params)

    @test nrow(df) > 0
    @test any(row -> row.name == :Vcmax25, eachrow(df))
    @test any(row -> row.name == :g1_w, eachrow(df))
  end
end
