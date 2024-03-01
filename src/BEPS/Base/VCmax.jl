function VCmax(lai, clumping, CosZs, param)
  # parameters for Vcmax-Nitrogen calculations
  G_theta = 0.5 # assuming a spherical leaf angle distribution
  K = G_theta * clumping / CosZs
  Kn = 0.3      # 0.713/2.4
  Vcmax0 = param[36+1]

  expr1 = 1 - exp(-K * lai)
  expr2 = 1 - exp(-lai * (Kn + K))
  expr3 = 1 - exp(-Kn * lai)

  # Formulas based on Chen et al., 2012, GBC
  if (expr1 > 0)
    Vcmax_sunlit = Vcmax0 * param[47+1] * param[46+1] * K * expr2 / (Kn + K) / expr1
  else
    Vcmax_sunlit = Vcmax0
  end

  if (K > 0 && lai > expr1 / K)
    Vcmax_shaded = Vcmax0 * param[47+1] * param[46+1] * (expr3 / Kn - expr2 / (Kn + K)) / (lai - expr1 / K)
  else
    Vcmax_shaded = Vcmax0
  end

  Vcmax_sunlit, Vcmax_shaded
end
