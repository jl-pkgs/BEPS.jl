"""
    readcoef(lc::Int, stxt::Int)

# Examples
```julia
lcs = [1, 6, 13, -1]
stxts = 1:11

coef1 = readcoef(1, 1)
```
"""
function readcoef(lc::Int, stxt::Int)
  coef = zeros(Cdouble, 48)
  
  clays = [0.03, 0.07, 0.1, 0.18, 0.15, 0.27, 0.34, 0.33, 0.4, 0.45, 0.6]
  clay_silts = [0.08, 0.19, 0.35, 0.58, 0.8, 0.4, 0.68, 0.91, 0.47, 0.9, 0.8]

  clay1 = clays[stxt]
  clay_silt1 = clay_silts[stxt]

  if lc in [1, 2, 3, 4, 5]  # conifer
    coef[1:8] = [0.301013, 0.148216, 0.212864, 0.347907, 0.02888, 0.02688, 0.1925, 0.5948]
    coef[32:48] = [650, 380, 100, 190, 350.26, 65.532, 56.532, 56.532, 56.532, 32.74, 350.26, 56.532, 56.532, 12, 10, 1.3557, 20]

    lignion_leaf = 0.302 / (0.15 + 0.018 * 0.6 * coef[37])
    lignion_fr = 0.280 / (0.15 + 0.018 * 0.6 * coef[38])
    lignion_wd = 0.40

    coef[9] = 3.9 * exp(-3 * lignion_leaf) * ((1 - lignion_leaf) * 0.6 + lignion_leaf * 0.3)
    coef[10] = 3.9 * exp(-3 * lignion_leaf) * (1 - lignion_leaf) * 0.4
    coef[11] = 3.9 * exp(-3 * lignion_leaf) * lignion_leaf * 0.7

    coef[12] = 14.8 * 0.6
    coef[13] = 14.8 * 0.4

    coef[14] = 4.8 * exp(-3 * lignion_fr) * ((1 - lignion_fr) * 0.55 + lignion_fr * 0.3)
    coef[15] = 4.8 * exp(-3 * lignion_fr) * (1 - lignion_fr) * 0.45
    coef[16] = 4.8 * exp(-3 * lignion_fr) * lignion_fr * 0.7

    coef[17] = 18.5 * 0.5
    coef[18] = 18.5 * 0.5

    coef[19] = 2.4 * exp(-3 * lignion_wd) * ((1 - lignion_wd) * 0.55 + lignion_wd * 0.45)
    coef[20] = 2.4 * exp(-3 * lignion_wd) * (1 - lignion_wd) * 0.45
    coef[21] = 2.4 * exp(-3 * lignion_wd) * lignion_wd * 0.55

    coef[22] = 7.3 * (1 - 0.75 * clay_silt1) * (0.85 - 0.68 * clay_silt1)
    coef[23] = 7.3 * (1 - 0.75 * clay_silt1) * (0.003 + 0.032 * clay1)
    coef[24] = 7.3 * (1 - 0.75 * clay_silt1) * (1 - (0.003 + 0.032 * clay1) - (0.85 - 0.68 * clay_silt1) - 5.0 / 18.0 * (0.01 + 0.04 * (1 - clay_silt1)))

    coef[25] = 6.0 * 0.6
    coef[26] = 6.0 * 0.4

    coef[27] = 0.25 * 0.55
    coef[28] = 0.25 * (0.003 - 0.009 * clay1)

    coef[28] = max(coef[28], 0.00001)
    if (0.003 - 0.009 * clay1) > 0.00001
      coef[29] = 0.25 * (1 - 0.55 - (0.003 - 0.009 * clay1))
    else
      coef[29] = 0.25 * (1 - 0.55 - 0.00001)
    end
    coef[30] = 0.007 * 0.5
    coef[31] = 0.007 * 0.5

  elseif lc in [6, 9]  # deciduous, tropic evergreen forest
    coef[1] = 0.422354
    coef[2] = 0.108994
    coef[3] = 0.242626
    coef[4] = 0.226026
    coef[5] = 0.01680
    coef[6] = 0.0248
    coef[7] = 1
    coef[8] = 0.5948

    coef[32] = 650
    coef[33] = 380
    coef[34] = 100
    coef[35] = 190
    coef[36] = 350.26
    coef[37] = 56.532
    coef[38] = 56.532
    coef[39] = 56.532
    coef[40] = 56.532
    coef[41] = 32.74
    coef[42] = 350.26
    coef[43] = 56.532
    coef[44] = 56.532
    coef[45] = 12.0
    coef[46] = 10.0
    coef[47] = 1.3557
    coef[48] = 20.0

    lignion_leaf = 0.224 / (0.15 + 0.018 * 0.6 * coef[37])
    lignion_fr = 0.200 / (0.15 + 0.018 * 0.6 * coef[38])
    lignion_wd = 0.30

    coef[9] = 3.9 * exp(-3 * lignion_leaf) * ((1 - lignion_leaf) * 0.6 + lignion_leaf * 0.3)
    coef[10] = 3.9 * exp(-3 * lignion_leaf) * (1 - lignion_leaf) * 0.4
    coef[11] = 3.9 * exp(-3 * lignion_leaf) * lignion_leaf * 0.7

    coef[12] = 14.8 * 0.6
    coef[13] = 14.8 * 0.4

    coef[14] = 4.8 * exp(-3 * lignion_fr) * ((1 - lignion_fr) * 0.55 + lignion_fr * 0.3)
    coef[15] = 4.8 * exp(-3 * lignion_fr) * (1 - lignion_fr) * 0.45
    coef[16] = 4.8 * exp(-3 * lignion_fr) * lignion_fr * 0.7

    coef[17] = 18.5 * 0.5
    coef[18] = 18.5 * 0.5

    coef[19] = 2.4 * exp(-3 * lignion_wd) * ((1 - lignion_wd) * 0.55 + lignion_wd * 0.45)
    coef[20] = 2.4 * exp(-3 * lignion_wd) * (1 - lignion_wd) * 0.45
    coef[21] = 2.4 * exp(-3 * lignion_wd) * lignion_wd * 0.55

    coef[22] = 7.3 * (1 - 0.75 * clay_silt1) * (0.85 - 0.68 * clay_silt1)
    coef[23] = 7.3 * (1 - 0.75 * clay_silt1) * (0.003 + 0.032 * clay1)
    coef[24] = 7.3 * (1 - 0.75 * clay_silt1) * (1 - (0.003 + 0.032 * clay1) - (0.85 - 0.68 * clay_silt1) - 5.0 / 18.0 * (0.01 + 0.04 * (1 - clay_silt1)))

    coef[25] = 6.0 * 0.6
    coef[26] = 6.0 * 0.4

    coef[27] = 0.25 * 0.55
    coef[28] = 0.25 * (0.003 - 0.009 * clay1)

    coef[28] = max(coef[28], 0.00001)
    if (0.003 - 0.009 * clay1) > 0.00001
      coef[29] = 0.25 * (1 - 0.55 - (0.003 - 0.009 * clay1))
    else
      coef[29] = 0.25 * (1 - 0.55 - 0.00001)
    end

    coef[30] = 0.0045 * 0.5
    coef[31] = 0.0045 * 0.5

  elseif lc == 13  # shrub

    coef[1] = 0.189428
    coef[2] = 0.053605
    coef[3] = 0.45
    coef[4] = 0.306967
    coef[5] = 0.025
    coef[6] = 0.04
    coef[7] = 0.8
    coef[8] = 0.75

    coef[32] = 650
    coef[33] = 380
    coef[34] = 100
    coef[35] = 190
    coef[36] = 370.26
    coef[37] = 63.532
    coef[38] = 63.532
    coef[39] = 63.532
    coef[40] = 63.532
    coef[41] = 32.74
    coef[42] = 370.26
    coef[43] = 63.532
    coef[44] = 63.532
    coef[45] = 12
    coef[46] = 10
    coef[47] = 1.3557
    coef[48] = 20

    lignion_leaf = 0.282 / (0.15 + 0.018 * 0.6 * coef[37])
    lignion_fr = 0.24 / (0.15 + 0.018 * 0.6 * coef[38])
    lignion_wd = 0.35

    coef[9] = 3.9 * exp(-3 * lignion_leaf) * ((1 - lignion_leaf) * 0.6 + lignion_leaf * 0.3)
    coef[10] = 3.9 * exp(-3 * lignion_leaf) * (1 - lignion_leaf) * 0.4
    coef[11] = 3.9 * exp(-3 * lignion_leaf) * lignion_leaf * 0.7

    coef[12] = 14.8 * 0.6
    coef[13] = 14.8 * 0.4

    coef[14] = 4.8 * exp(-3 * lignion_fr) * ((1 - lignion_fr) * 0.55 + lignion_fr * 0.3)
    coef[15] = 4.8 * exp(-3 * lignion_fr) * (1 - lignion_fr) * 0.45
    coef[16] = 4.8 * exp(-3 * lignion_fr) * lignion_fr * 0.7

    coef[17] = 18.5 * 0.5
    coef[18] = 18.5 * 0.5

    coef[19] = 2.4 * exp(-3 * lignion_wd) * ((1 - lignion_wd) * 0.55 + lignion_wd * 0.45)
    coef[20] = 2.4 * exp(-3 * lignion_wd) * (1 - lignion_wd) * 0.45
    coef[21] = 2.4 * exp(-3 * lignion_wd) * lignion_wd * 0.55

    coef[22] = 7.3 * (1 - 0.75 * clay_silt1) * (0.85 - 0.68 * clay_silt1)
    coef[23] = 7.3 * (1 - 0.75 * clay_silt1) * (0.003 + 0.032 * clay1)
    coef[24] = 7.3 * (1 - 0.75 * clay_silt1) * (1 - (0.003 + 0.032 * clay1) - (0.85 - 0.68 * clay_silt1) - 5.0 / 18.0 * (0.01 + 0.04 * (1 - clay_silt1)))

    coef[25] = 6.0 * 0.6
    coef[26] = 6.0 * 0.4

    coef[27] = 0.25 * 0.55
    coef[28] = 0.25 * (0.003 - 0.009 * clay1)

    coef[28] = max(coef[28], 0.00001)
    if (0.003 - 0.009 * clay1) > 0.00001
      coef[29] = 0.25 * (1 - 0.55 - (0.003 - 0.009 * clay1))
    else
      coef[29] = 0.25 * (1 - 0.55 - 0.00001)
    end

    coef[30] = 0.007 * 0.5
    coef[31] = 0.007 * 0.5
  else

    coef[1] = 0.331684
    coef[2] = 0.053605
    coef[3] = 0.307745
    coef[4] = 0.306967
    coef[5] = 0.0278
    coef[6] = 0.0448
    coef[7] = 0.39448
    coef[8] = 0.5948

    coef[32] = 650
    coef[33] = 380
    coef[34] = 100
    coef[35] = 190
    coef[36] = 370.26
    coef[37] = 63.532
    coef[38] = 63.532
    coef[39] = 63.532
    coef[40] = 63.532
    coef[41] = 32.74
    coef[42] = 370.26
    coef[43] = 63.532
    coef[44] = 63.532
    coef[45] = 12
    coef[46] = 10
    coef[47] = 1.3557
    coef[48] = 20

    lignion_leaf = 0.6 * 0.224 / (0.15 + 0.018 * 0.6 * coef[37])
    lignion_fr = 0.6 * 0.200 / (0.15 + 0.018 * 0.6 * coef[38])
    lignion_wd = 0.30

    coef[9] = 3.9 * exp(-3 * lignion_leaf) * ((1 - lignion_leaf) * 0.6 + lignion_leaf * 0.3)
    coef[10] = 3.9 * exp(-3 * lignion_leaf) * (1 - lignion_leaf) * 0.4
    coef[11] = 3.9 * exp(-3 * lignion_leaf) * lignion_leaf * 0.7

    coef[12] = 14.8 * 0.6
    coef[13] = 14.8 * 0.4

    coef[14] = 4.8 * exp(-3 * lignion_fr) * ((1 - lignion_fr) * 0.55 + lignion_fr * 0.3)
    coef[15] = 4.8 * exp(-3 * lignion_fr) * (1 - lignion_fr) * 0.45
    coef[16] = 4.8 * exp(-3 * lignion_fr) * lignion_fr * 0.7

    coef[17] = 18.5 * 0.5
    coef[18] = 18.5 * 0.5

    coef[19] = 2.4 * exp(-3 * lignion_wd) * ((1 - lignion_wd) * 0.55 + lignion_wd * 0.45)
    coef[20] = 2.4 * exp(-3 * lignion_wd) * (1 - lignion_wd) * 0.45
    coef[21] = 2.4 * exp(-3 * lignion_wd) * lignion_wd * 0.55

    coef[22] = 7.3 * (1 - 0.75 * clay_silt1) * (0.85 - 0.68 * clay_silt1)
    coef[23] = 7.3 * (1 - 0.75 * clay_silt1) * (0.003 + 0.032 * clay1)
    coef[24] = 7.3 * (1 - 0.75 * clay_silt1) * (1 - (0.003 + 0.032 * clay1) - (0.85 - 0.68 * clay_silt1) - 5.0 / 18.0 * (0.01 + 0.04 * (1 - clay_silt1)))

    coef[25] = 6.0 * 0.6
    coef[26] = 6.0 * 0.4

    coef[27] = 0.25 * 0.55
    coef[28] = 0.25 * (0.003 - 0.009 * clay1)

    coef[28] = max(coef[28], 0.00001)
    if (0.003 - 0.009 * clay1) > 0.00001
      coef[29] = 0.25 * (1 - 0.55 - (0.003 - 0.009 * clay1))
    else
      coef[29] = 0.25 * (1 - 0.55 - 0.00001)
    end

    coef[30] = 0.007 * 0.5
    coef[31] = 0.007 * 0.5
  end

  coef
end
