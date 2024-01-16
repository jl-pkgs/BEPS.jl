using Test
using BEPS
# import BEPS: readparam, readcoef

@testset "readparam" begin
  lcs = [1, 2, 6, 9, 13, 40, -1]

  for lc = lcs
    par1 = readparam(lc)
    par2 = clang.readparam(lc)
    @test maximum(abs.(par1 .- par2)) == 0
  end
end


@testset "readparam" begin
  lcs = [1, 6, 13, -1]
  stxts = 1:11

  for lc = lcs, stxt = stxts
    coef1 = readcoef(lc, stxt)
    coef2 = clang.readcoef(lc, stxt)
    @test maximum(abs.(coef1 .- coef2)) == 0
  end
end
