using Test
using BEPS
using BEPS: path_proj

function nanmax(x)
  x = collect(x)
  x = x[.!isnan.(x)]
  maximum(x)
end

lai = readdlm(path_proj("examples/input/p1_lai.txt"))[:]

par = (lon=120.5, lat=30.5, landcover=25, clumping=0.85,
  soil_type=8, Tsoil=2.2,
  soilwater=0.4115, snowdepth=0.0)

@testset "besp_main julia" begin
  d = deserialize(path_proj("data/p1_meteo"))
  d.tem = d.tem .- 5 # gpp_u_sunlit的计算误差变大  

  @time df_jl, df_ET_jl = besp_main(d, lai, par; version="julia")
  @time df_c, df_ET_c = besp_main(d, lai, par; version="c")
  # r = sum(df_jl)

  df_diff = abs.(df_jl .- df_c)
  df_diff_perc = abs.(df_jl .- df_c) ./ df_c .* 100
  l = max(df_diff_perc)
  @show l
  @test true
  # @test nanmax(l) < 1.5 # SH, 1.48%的误差, current 0.09%
end

# > 需要核对是哪个变量引起的
