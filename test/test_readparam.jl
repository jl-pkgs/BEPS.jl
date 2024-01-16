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


# Init_Soil_Parameters(par.landcover, par.soil_type, parameter[28], p_soil)
@testset "Init_Soil_Parameters" begin
  p_jl = Soil()
  p_c = Soil_c()

  for stxt = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    Init_Soil_Parameters(1, stxt, 0.1, p_jl)
    clang.Init_Soil_Parameters(1, stxt, 0.1, p_c)

    names = fieldnames(typeof(p_jl))
    # name = names[1]
    for name in names
      x_jl = getfield(p_jl, name)
      x_c = getfield(p_c, name)

      # println(name)
      @test maximum(abs.(x_c .- x_jl)) <= 1e-8
    end
  end
end

# name = :r_root_decay
