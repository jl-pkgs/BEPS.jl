using Test


@testset "readparam" begin
  lcs = [1, 2, 6, 9, 13, 40, -1]

  for lc = lcs
    par1 = readparam(lc)
    par2 = clang.readparam(lc)
    @test maximum(abs.(par1 .- par2)) == 0
  end
end


