# VCmax-Nitrogen calculations
# 
# VCmax25 = 62.5          # maximum capacity of Rubisco at 25C-VCmax	
# N_leaf = 3.10 + 1.35    # leaf Nitrogen content	mean value + 1 SD g/m2 
# slope = 20.72 / 62.5       # slope of VCmax-N curve
function VCmax(lai::FT, Ω::FT, CosZs::FT, VCmax25::FT, N_leaf::FT, χ::FT) where {FT<:Real}
  CosZs <= 0 && return 0.0, 0.0 # 光合仅发生在白天
  K = 0.5 * Ω / CosZs # assuming a spherical leaf angle distribution
  Kn = 0.3            # 0.713/2.4

  expr1 = 1 - exp(-K * lai)
  expr2 = 1 - exp(-lai * (Kn + K))
  expr3 = 1 - exp(-Kn * lai)

  # Formulas based on Chen et al., 2012, GBC
  if (expr1 > 0)
    VCmax_sunlit = VCmax25 * χ * N_leaf * K * expr2 / (Kn + K) / expr1
  else
    VCmax_sunlit = VCmax25
  end

  if (K > 0 && lai > expr1 / K)
    VCmax_shaded = VCmax25 * χ * N_leaf *
                   (expr3 / Kn - expr2 / (Kn + K)) / (lai - expr1 / K)
  else
    VCmax_shaded = VCmax25
  end
  VCmax_sunlit, VCmax_shaded
end
