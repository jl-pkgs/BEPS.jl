function lai2(clumping::Float64, CosZs::Float64,
  stem_o::Float64, stem_u::Float64,
  lai_o::Float64, lai_u::Float64,
  LAI::Leaf, PAI::Leaf)

  PAI.o_sunlit = CosZs > 0 ? 2 * CosZs * (1 - exp(-0.5 * clumping * (lai_o + stem_o) / CosZs)) : 0
  PAI.o_shaded = (lai_o + stem_o) - PAI.o_sunlit

  PAI.u_sunlit = CosZs > 0 ? 2 * CosZs * (1 - exp(-0.5 * clumping * (lai_o + stem_o + lai_u + stem_u) / CosZs)) - PAI.o_sunlit : 0
  PAI.u_shaded = (lai_u + stem_u) - PAI.u_sunlit

  LAI.o_sunlit = CosZs > 0 ? 2 * CosZs * (1 - exp(-0.5 * clumping * lai_o / CosZs)) : 0
  LAI.o_shaded = max(0, lai_o - LAI.o_sunlit)  # edited by J. Leng

  LAI.u_sunlit = CosZs > 0 ? 2 * CosZs * (1 - exp(-0.5 * clumping * (lai_o + lai_u) / CosZs)) - LAI.o_sunlit : 0
  LAI.u_shaded = max(0, lai_u - LAI.u_sunlit)  # edited by J. Leng
end

function lai2(clumping::Float64, CosZs::Float64,
  stem_o::Float64, stem_u::Float64,
  lai_o::Float64, lai_u::Float64)

  LAI = Leaf()
  PAI = Leaf()
  
  lai2(clumping, CosZs, stem_o, stem_u, lai_o, lai_u, LAI, PAI)
  LAI, PAI
end
