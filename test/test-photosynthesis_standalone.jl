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
  end
  
  @testset "基本功能" begin
    params = InitParam_Photo_Farquhar()
    
    result = BEPS.Photosynthesis.photosynthesis(
      25.0, 60.0, 500.0, 2.0, params
    )
    
    @test result.An.sunlit > 0
    @test result.An.shaded >= 0
    @test result.Gs.sunlit > 0
    @test result.Gs.shaded >= 0
    @test result.An.sunlit > result.An.shaded
  end
  
  @testset "温度响应" begin
    params = InitParam_Photo_Farquhar()
    
    result_cold = BEPS.Photosynthesis.photosynthesis(10.0, 60.0, 500.0, 2.0, params)
    result_warm = BEPS.Photosynthesis.photosynthesis(25.0, 60.0, 500.0, 2.0, params)
    result_hot = BEPS.Photosynthesis.photosynthesis(40.0, 60.0, 500.0, 2.0, params)
    
    @test result_warm.An.sunlit > result_cold.An.sunlit
    @test result_hot.An.sunlit < result_warm.An.sunlit
  end
  
  @testset "光响应" begin
    params = InitParam_Photo_Farquhar()
    
    result_dark = BEPS.Photosynthesis.photosynthesis(25.0, 60.0, 0.0, 2.0, params)
    result_light = BEPS.Photosynthesis.photosynthesis(25.0, 60.0, 500.0, 2.0, params)
    
    @test result_light.An.sunlit > result_dark.An.sunlit
    @test result_dark.An.sunlit < 0  # 暗呼吸占主导
  end
  
  @testset "水分胁迫" begin
    params = InitParam_Photo_Farquhar()
    
    result_no_stress = BEPS.Photosynthesis.photosynthesis(25.0, 60.0, 500.0, 2.0, params; β_soil=1.0)
    result_stress = BEPS.Photosynthesis.photosynthesis(25.0, 60.0, 500.0, 2.0, params; β_soil=0.3)
    
    @test result_no_stress.An.sunlit > result_stress.An.sunlit
  end
  
  @testset "参数堆叠" begin
    params = InitParam_Photo_Farquhar()
    df = parameters(params)
    
    @test nrow(df) > 0
    @test any(row -> row.name == :Vcmax25, eachrow(df))
    @test any(row -> row.name == :g1_w, eachrow(df))
  end
end
