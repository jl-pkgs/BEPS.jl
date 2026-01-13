using Test
using BEPS

# 误差主要集中: slope_Vc
@testset "readVegParam" begin
  lcs = [1, 2, 6, 9, 13, 40]
  for lc = lcs
    par1 = readVegParam(lc)
    par2 = clang.readVegParam(lc)
    # 对于-1, julia修改了LAI_max_understory，from [0.01] to [2.4]
    @test par1 ≈ par2 atol = 1e-3
  end
end
# lc = -1

# @testset "readcoef" begin
#   lcs = [1, 6, 13, -1]
#   stxts = 1:11
#   for lc = lcs, stxt = stxts
#     coef1 = readcoef(lc, stxt)
#     coef2 = clang.readcoef(lc, stxt)
#     @test maximum(abs.(coef1 .- coef2)) == 0
#   end
# end
