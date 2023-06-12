function netRadiation(shortRad_global::Float64, CosZs::Float64, temp_o::Float64, temp_u::Float64, temp_g::Float64,
  lai_o::Float64, lai_u::Float64, lai_os::Float64, lai_us::Float64,
  lai::Float64,
  clumping::Float64, temp_air::Float64, rh::Float64,
  albedo_snow_v::Float64, albedo_snow_n::Float64, percentArea_snow_o::Float64, percentArea_snow_u::Float64, percent_snow_g::Float64,
  albedo_v_o::Float64, albedo_n_o::Float64, albedo_v_u::Float64, albedo_n_u::Float64, albedo_v_g::Float64, albedo_n_g::Float64,
  # netRad_o::Float64, netRad_u::Float64, netRad_g::Float64,
  netRadLeaf::Leaf,
  netShortRadLeaf::Leaf)

  # 计算剩余代码...
  sb_constant = 5.67 / 100000000                            # stephen-boltzman constant

  # calculate albedo of canopy in this step
  albedo_v_os = albedo_v_o * (1 - percentArea_snow_o) + albedo_snow_v * percentArea_snow_o  # visible, overstory
  albedo_n_os = albedo_n_o * (1 - percentArea_snow_o) + albedo_snow_n * percentArea_snow_o  # near infrared
  albedo_v_us = albedo_v_u * (1 - percentArea_snow_u) + albedo_snow_v * percentArea_snow_u  #        , understory
  albedo_n_us = albedo_n_u * (1 - percentArea_snow_u) + albedo_snow_n * percentArea_snow_u

  albedo_o = 0.5 * (albedo_v_os + albedo_n_os)
  albedo_u = 0.5 * (albedo_v_us + albedo_n_us)

  # calculate albedo of ground in this step
  albedo_v_gs = albedo_v_g * (1 - percent_snow_g) + albedo_snow_v * percent_snow_g
  albedo_n_gs = albedo_n_g * (1 - percent_snow_g) + albedo_snow_n * percent_snow_g
  albedo_g = 0.5 * (albedo_v_gs + albedo_n_gs)

  # separate global solar radiation into direct and diffuse one
  if (CosZs < 0.001)  # solar zenith angle small, all diffuse radiation
    ratio_cloud = 0
  else
    ratio_cloud = shortRad_global / (1367 * CosZs)  # Luo2018, A4
  end

  if (ratio_cloud > 0.8)
    shortRad_df = 0.13 * shortRad_global  # Luo2018, A2
  else
    shortRad_df = (0.943 + 0.734 * ratio_cloud - 4.9 * pow((ratio_cloud), 2) +
                   1.796 * pow((ratio_cloud), 3) + 2.058 * pow((ratio_cloud), 4)) * shortRad_global  # Luo2018, A2
  end
  shortRad_df = clamp(shortRad_df, 0, shortRad_global)
  shortRad_dir = shortRad_global - shortRad_df  # Luo2018, A3

  # fraction at each layer of canopy, direct and diffuse. use Leaf only lai here
  gap_o_dir = exp(-0.5 * clumping * lai_o / CosZs)
  gap_u_dir = exp(-0.5 * clumping * lai_u / CosZs)

  # indicators to describe leaf distribution angles in canopy. slightly related with LAI
  cosQ_o = 0.537 + 0.025 * lai_o  # Luo2018, A10, a representative zenith angle for diffuse radiation transmission
  cosQ_u = 0.537 + 0.025 * lai_u

  gap_o_df = exp(-0.5 * clumping * lai_o / cosQ_o)
  gap_u_df = exp(-0.5 * clumping * lai_u / cosQ_u)

  gap_os_df = exp(-0.5 * clumping * lai_os / cosQ_o)  # considering stem
  gap_us_df = exp(-0.5 * clumping * lai_us / cosQ_u)

  # emissivity of each part
  e_actual = cal_ea(temp_air, rh)

  emissivity_air = 1 - exp(-(pow(e_actual * 10.0, (temp_air + 273.15) / 1200.0)))
  emissivity_air = cal_ea(emissivity_air, 0.7, 1.0)

  emissivity_o = 0.98
  emissivity_u = 0.98
  emissivity_g = 0.96

  if shortRad_global > zero && CosZs > zero
    # net short direct radiation on canopy and ground
    netShortRad_o_dir = shortRad_dir * ((1 - albedo_o) - (1 - albedo_u) * gap_o_dir)  # dir into dif_under
    netShortRad_u_dir = shortRad_dir * gap_o_dir * ((1 - albedo_u) - (1 - albedo_g) * gap_u_dir)
    netShortRad_g_dir = shortRad_dir * gap_o_dir * gap_u_dir * (1 - albedo_g)
  else
    netShortRad_o_dir = 0
    netShortRad_u_dir = 0
    netShortRad_g_dir = 0
  end

  if shortRad_global > zero && CosZs > zero
    # net short diffuse radiation on canopy and ground
    netShortRad_o_df = shortRad_df * ((1 - albedo_o) - (1 - albedo_u) * gap_o_df) +
                       0.21 * clumping * shortRad_dir * (1.1 - 0.1 * lai_o) * exp(-CosZs)  # A8
    netShortRad_u_df = shortRad_df * gap_o_df * ((1 - albedo_u) - (1 - albedo_g) * gap_u_df) +
                       0.21 * clumping * shortRad_dir * gap_o_dir * (1.1 - 0.1 * lai_u) * exp(-CosZs)  # A9
    netShortRad_g_df = shortRad_df * gap_o_df * gap_u_df * (1 - albedo_g)
  else
    netShortRad_o_df = 0
    netShortRad_u_df = 0
    netShortRad_g_df = 0
  end

  # total net shortwave radiation at canopy level
  netShortRad_o = netShortRad_o_dir + netShortRad_o_df
  netShortRad_u = netShortRad_u_dir + netShortRad_u_df
  netShortRad_g = netShortRad_g_dir + netShortRad_g_df

  # 计算植被和地面的净长波辐射
  longRad_air = emissivity_air * sb_constant * (temp_air + 273.15)^4
  longRad_o = emissivity_o * sb_constant * (temp_o + 273.15)^4
  longRad_u = emissivity_u * sb_constant * (temp_u + 273.15)^4
  longRad_g = emissivity_g * sb_constant * (temp_g + 273.15)^4

  netLongRad_o = (emissivity_o * (longRad_air + longRad_u * (1 - gap_u_df) + longRad_g * gap_u_df) - 2 * longRad_o) *
                 (1 - gap_o_df) +
                 emissivity_o * (1 - emissivity_u) * (1 - gap_u_df) * (longRad_air * gap_o_df + longRad_o * (1 - gap_o_df))

  netLongRad_u = (emissivity_u * (longRad_air * gap_o_df + longRad_o * (1 - gap_o_df) + longRad_g) - 2 * longRad_u) * (1 - gap_u_df) +
                 (1 - emissivity_g) * ((longRad_air * gap_o_df + longRad_o * (1 - gap_o_df)) * gap_u_df + longRad_u * (1 - gap_u_df)) +
                 emissivity_u * (1 - emissivity_o) * (longRad_u * (1 - gap_u_df) + longRad_g * gap_u_df) * (1 - gap_o_df)

  netLongRad_g = emissivity_g * ((longRad_air * gap_o_df + longRad_o * (1 - gap_o_df)) * gap_u_df + longRad_u * (1 - gap_u_df)) -
                 longRad_g + (1 - emissivity_u) * longRad_g * (1 - gap_u_df)

  # 计算植被和地面的总净辐射
  netRad_o = netShortRad_o + netLongRad_o
  netRad_u = netShortRad_u + netLongRad_u
  netRad_g = netShortRad_g + netLongRad_g


  if shortRad_global > zero && CosZs > zero # only happens in day time, when sun is out
    shortRadLeaf_o_dir = 0.5 * shortRad_dir / CosZs
    shortRadLeaf_o_dir = min(shortRadLeaf_o_dir, 0.7 * 1362)
    shortRadLeaf_u_dir = shortRadLeaf_o_dir

    shortRadLeaf_o_df = (shortRad_df - shortRad_df * gap_os_df) / lai_os + 0.07 * shortRad_dir * (1.1 - 0.1 * lai_os) * exp(-CosZs)
    shortRadLeaf_u_df = (shortRad_df * gap_o_df - shortRad_df * gap_o_df * gap_us_df) / lai_us + 0.05 * shortRad_dir * gap_o_dir * (1.1 - 0.1 * lai_us) * exp(-CosZs)
  else
    shortRadLeaf_o_dir = 0
    shortRadLeaf_u_dir = 0
    shortRadLeaf_o_df = 0
    shortRadLeaf_u_df = 0
  end

  # overstorey sunlit leaves
  if lai.o_sunlit > 0
    netShortRadLeaf.o_sunlit = (shortRadLeaf_o_dir + shortRadLeaf_o_df) * (1 - albedo_o) # diffuse
    netLongRadLeaf.o_sunlit = netLongRad_o / lai_os # leaf level net long
    netRadLeaf.o_sunlit = netShortRadLeaf.o_sunlit + netLongRadLeaf.o_sunlit
  else
    netShortRadLeaf.o_sunlit = (shortRadLeaf_o_dir + shortRadLeaf_o_df) * (1 - albedo_o)
    netLongRadLeaf.o_sunlit = netLongRad_o
    netRadLeaf.o_sunlit = netShortRadLeaf.o_sunlit + netLongRadLeaf.o_sunlit
  end

  # overstorey shaded leaf
  if lai.o_shaded > 0
    netShortRadLeaf.o_shaded = shortRadLeaf_o_df * (1 - albedo_o) # diffuse
    netLongRadLeaf.o_shaded = netLongRad_o / lai_os

    netRadLeaf.o_shaded = netShortRadLeaf.o_shaded + netLongRadLeaf.o_shaded
  else
    netShortRadLeaf.o_shaded = shortRadLeaf_o_df * (1 - albedo_o) # diffuse
    netLongRadLeaf.o_shaded = netLongRad_o
    netRadLeaf.o_shaded = netShortRadLeaf.o_shaded + netLongRadLeaf.o_shaded
  end

  # understorey sunlit leaf
  if lai.u_sunlit > 0
    netShortRadLeaf.u_sunlit = (shortRadLeaf_u_dir + shortRadLeaf_u_df) * (1 - albedo_u)
    netLongRadLeaf.u_sunlit = netLongRad_u / lai_us

    netRadLeaf.u_sunlit = netShortRadLeaf.u_sunlit + netLongRadLeaf.u_sunlit
  else
    netShortRadLeaf.u_sunlit = (shortRadLeaf_u_dir + shortRadLeaf_u_df) * (1 - albedo_u)
    netLongRadLeaf.u_sunlit = netLongRad_u

    netRadLeaf.u_sunlit = netShortRadLeaf.u_sunlit + netLongRadLeaf.u_sunlit
  end

  if lai.u_shaded > 0
    netShortRadLeaf.u_shaded = shortRadLeaf_u_df * (1 - albedo_u)
    netLongRadLeaf.u_shaded = netLongRad_u / lai_us

    netRadLeaf.u_shaded = netShortRadLeaf.u_shaded + netLongRadLeaf.u_shaded
  else
    netShortRadLeaf.u_shaded = shortRadLeaf_u_df * (1 - albedo_u)
    netLongRadLeaf.u_shaded = netLongRad_u

    netRadLeaf.u_shaded = netShortRadLeaf.u_shaded + netLongRadLeaf.u_shaded
  end

  # 叶片尺度的净辐射更新方式
  # 参考Chen 2012年的聚集指数论文
  return netRadLeaf
end

