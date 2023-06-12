function netRadiation_c(shortRad_global, CosZs, 
  temp_o, temp_u, temp_g,
  lai_o, lai_u, lai_os, lai_us, lai::Leaf, clumping, temp_air, rh,
  albedo_snow_v, albedo_snow_n,
  percentArea_snow_o, percentArea_snow_u, percent_snow_g,
  albedo_v_o, albedo_n_o, albedo_v_u, albedo_n_u, albedo_v_g, albedo_n_g,
  # netRad_o::TypeRef, netRad_u::TypeRef, netRad_g::TypeRef,
  netRadLeaf::Leaf, netShortRadLeaf::Leaf, Rnl_Leaf::Leaf)

  netRad_o = init_dbl()
  netRad_u = init_dbl()
  netRad_g = init_dbl()

  ccall((:netRadiation, libbeps), Cvoid,
    (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Leaf, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble,
      Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Leaf}, Ptr{Leaf}),
    shortRad_global, CosZs, temp_o, temp_u, temp_g, lai_o, lai_u, lai_os, lai_us, lai, clumping, temp_air, rh, albedo_snow_v, albedo_snow_n, percentArea_snow_o, percentArea_snow_u, percent_snow_g,
    albedo_v_o, albedo_n_o, albedo_v_u, albedo_n_u, albedo_v_g, albedo_n_g,
    netRad_o, netRad_u, netRad_g, Ref(netRadLeaf), Ref(netShortRadLeaf))
  
  netRad_o[], netRad_u[], netRad_g[]
end


function netRadiation_jl(Rs_global::T, CosZs::T,
  temp_o::T, temp_u::T, temp_g::T,
  lai_o::T, lai_u::T, lai_os::T, lai_us::T,
  lai::Leaf,
  clumping::T, temp_air::T, rh::T,
  albedo_snow_v::T, albedo_snow_n::T,
  percArea_snow_o::T, percArea_snow_u::T, perc_snow_g::T,
  albedo_v_o::T, albedo_n_o::T, albedo_v_u::T, albedo_n_u::T, 
  albedo_v_g::T, albedo_n_g::T,
  # Rn_o::T, Rn_u::T, Rn_g::T,
  Rn_Leaf::Leaf, 
  Rns_Leaf::Leaf, 
  Rnl_Leaf::Leaf) where {T<:Real}
  # Rnl_Leaf::Leaf = Leaf()
  
  # calculate albedo of canopy in this step
  albedo_v_os = albedo_v_o * (1 - percArea_snow_o) + albedo_snow_v * percArea_snow_o  # visible, overstory
  albedo_n_os = albedo_n_o * (1 - percArea_snow_o) + albedo_snow_n * percArea_snow_o  # near infrared
  albedo_v_us = albedo_v_u * (1 - percArea_snow_u) + albedo_snow_v * percArea_snow_u  # , understory
  albedo_n_us = albedo_n_u * (1 - percArea_snow_u) + albedo_snow_n * percArea_snow_u

  albedo_o = 0.5 * (albedo_v_os + albedo_n_os)
  albedo_u = 0.5 * (albedo_v_us + albedo_n_us)

  # calculate albedo of ground in this step
  albedo_v_gs = albedo_v_g * (1 - perc_snow_g) + albedo_snow_v * perc_snow_g
  albedo_n_gs = albedo_n_g * (1 - perc_snow_g) + albedo_snow_n * perc_snow_g
  albedo_g = 0.5 * (albedo_v_gs + albedo_n_gs)

  # separate global solar radiation into direct and diffuse one
  if (CosZs < 0.001)  # solar zenith angle small, all diffuse radiation
    ratio_cloud = 0.0
  else
    ratio_cloud = Rs_global / (1367 * CosZs)  # Luo2018, A4
  end

  if (ratio_cloud > 0.8)
    Rs_df = 0.13 * Rs_global  # Luo2018, A2
  else
    Rs_df = (0.943 + 0.734 * ratio_cloud - 4.9 * pow((ratio_cloud), 2) +
             1.796 * pow((ratio_cloud), 3) + 2.058 * pow((ratio_cloud), 4)) * Rs_global  # Luo2018, A2
  end

  Rs_df = clamp(Rs_df, 0.0, Rs_global)
  Rs_dir = Rs_global - Rs_df  # Luo2018, A3

  # fraction at each layer of canopy, direct and diffuse. use Leaf only lai here
  gap_o_dir = exp(-0.5 * clumping * lai_o / CosZs)
  gap_u_dir = exp(-0.5 * clumping * lai_u / CosZs)

  # indicators to describe leaf distribution angles in canopy. slightly related with LAI
  cosQ_o = 0.537 + 0.025 * lai_o  # Luo2018, A10, a representative zenith angle for diffuse radiation transmission
  cosQ_u = 0.537 + 0.025 * lai_u

  gap_o_df = exp(-0.5 * clumping * lai_o / cosQ_o)
  gap_os_df = exp(-0.5 * clumping * lai_os / cosQ_o)  # considering stem

  gap_u_df = exp(-0.5 * clumping * lai_u / cosQ_u)
  gap_us_df = exp(-0.5 * clumping * lai_us / cosQ_u)

  # emissivity of each part
  ea = cal_ea(temp_air, rh)
  emiss_air = 1 - exp(-(pow(ea * 10.0, (temp_air + 273.15) / 1200.0)))
  emiss_air = clamp(emiss_air, 0.7, 1.0)

  emiss_o = 0.98
  emiss_u = 0.98
  emiss_g = 0.96

  if Rs_global > zero && CosZs > zero
    # net short direct radiation on canopy and ground
    Rns_o_dir = Rs_dir * ((1 - albedo_o) - (1 - albedo_u) * gap_o_dir)  # dir into dif_under
    Rns_u_dir = Rs_dir * gap_o_dir * ((1 - albedo_u) - (1 - albedo_g) * gap_u_dir)
    Rns_g_dir = Rs_dir * gap_o_dir * gap_u_dir * (1 - albedo_g)
  else
    Rns_o_dir = 0.0
    Rns_u_dir = 0.0
    Rns_g_dir = 0.0
  end

  if Rs_global > zero && CosZs > zero
    # net short diffuse radiation on canopy and ground
    Rns_o_df = Rs_df * ((1 - albedo_o) - (1 - albedo_u) * gap_o_df) +
               0.21 * clumping * Rs_dir * (1.1 - 0.1 * lai_o) * exp(-CosZs)  # A8
    Rns_u_df = Rs_df * gap_o_df * ((1 - albedo_u) - (1 - albedo_g) * gap_u_df) +
               0.21 * clumping * Rs_dir * gap_o_dir * (1.1 - 0.1 * lai_u) * exp(-CosZs)  # A9
    Rns_g_df = Rs_df * gap_o_df * gap_u_df * (1 - albedo_g)
  else
    Rns_o_df = 0.0
    Rns_u_df = 0.0
    Rns_g_df = 0.0
  end

  # total net shortwave radiation at canopy level
  Rns_o = Rns_o_dir + Rns_o_df
  Rns_u = Rns_u_dir + Rns_u_df
  Rns_g = Rns_g_dir + Rns_g_df

  # 计算植被和地面的净长波辐射
  Rl_air = cal_Rln_out(emiss_air, temp_air)
  Rl_o = cal_Rln_out(emiss_o, temp_o)
  Rl_u = cal_Rln_out(emiss_u, temp_u)
  Rl_g = cal_Rln_out(emiss_g, temp_g)

  Rnl_o = (emiss_o * (Rl_air + Rl_u * (1 - gap_u_df) + Rl_g * gap_u_df) - 2 * Rl_o) *
          (1 - gap_o_df) +
          emiss_o * (1 - emiss_u) * (1 - gap_u_df) * (Rl_air * gap_o_df + Rl_o * (1 - gap_o_df))

  Rnl_u = (emiss_u * (Rl_air * gap_o_df + Rl_o * (1 - gap_o_df) + Rl_g) - 2 * Rl_u) * (1 - gap_u_df) +
          (1 - emiss_g) * ((Rl_air * gap_o_df + Rl_o * (1 - gap_o_df)) * gap_u_df + Rl_u * (1 - gap_u_df)) +
          emiss_u * (1 - emiss_o) * (Rl_u * (1 - gap_u_df) + Rl_g * gap_u_df) * (1 - gap_o_df)

  Rnl_g = emiss_g * ((Rl_air * gap_o_df + Rl_o * (1 - gap_o_df)) * gap_u_df + Rl_u * (1 - gap_u_df)) -
          Rl_g + (1 - emiss_u) * Rl_g * (1 - gap_u_df)

  # 计算植被和地面的总净辐射
  Rn_o = Rns_o + Rnl_o
  Rn_u = Rns_u + Rnl_u
  Rn_g = Rns_g + Rnl_g

  if Rs_global > zero && CosZs > zero # only happens in day time, when sun is out
    Rs_o_dir = 0.5 * Rs_dir / CosZs
    Rs_o_dir = min(Rs_o_dir, 0.7 * 1362)
    Rs_u_dir = Rs_o_dir

    Rs_o_df = (Rs_df - Rs_df * gap_os_df) / lai_os + 0.07 * Rs_dir * (1.1 - 0.1 * lai_os) * exp(-CosZs)
    Rs_u_df = (Rs_df * gap_o_df - Rs_df * gap_o_df * gap_us_df) / lai_us +
                  0.05 * Rs_dir * gap_o_dir * (1.1 - 0.1 * lai_us) * exp(-CosZs)
  else
    Rs_o_dir = 0.
    Rs_u_dir = 0.
    Rs_o_df = 0.
    Rs_u_df = 0.
  end

  # overstorey sunlit leaves
  if lai.o_sunlit > 0
    Rns_Leaf.o_sunlit = (Rs_o_dir + Rs_o_df) * (1 - albedo_o) # diffuse
    Rnl_Leaf.o_sunlit = Rnl_o / lai_os # leaf level net long
  else
    Rns_Leaf.o_sunlit = (Rs_o_dir + Rs_o_df) * (1 - albedo_o)
    Rnl_Leaf.o_sunlit = Rnl_o
  end

  # overstorey shaded leaf
  if lai.o_shaded > 0
    Rns_Leaf.o_shaded = Rs_o_df * (1 - albedo_o) # diffuse
    Rnl_Leaf.o_shaded = Rnl_o / lai_os
  else
    Rns_Leaf.o_shaded = Rs_o_df * (1 - albedo_o) # diffuse
    Rnl_Leaf.o_shaded = Rnl_o
  end
  
  # understorey sunlit leaf
  if lai.u_sunlit > 0
    Rns_Leaf.u_sunlit = (Rs_u_dir + Rs_u_df) * (1 - albedo_u)
    Rnl_Leaf.u_sunlit = Rnl_u / lai_us
  else
    Rns_Leaf.u_sunlit = (Rs_u_dir + Rs_u_df) * (1 - albedo_u)
    Rnl_Leaf.u_sunlit = Rnl_u
  end
  
  if lai.u_shaded > 0
    Rns_Leaf.u_shaded = Rs_u_df * (1 - albedo_u)
    Rnl_Leaf.u_shaded = Rnl_u / lai_us
  else
    Rns_Leaf.u_shaded = Rs_u_df * (1 - albedo_u)
    Rnl_Leaf.u_shaded = Rnl_u
  end

  Rn_Leaf.o_sunlit = Rns_Leaf.o_sunlit + Rnl_Leaf.o_sunlit
  Rn_Leaf.o_shaded = Rns_Leaf.o_shaded + Rnl_Leaf.o_shaded
  Rn_Leaf.u_sunlit = Rns_Leaf.u_sunlit + Rnl_Leaf.u_sunlit
  Rn_Leaf.u_shaded = Rns_Leaf.u_shaded + Rnl_Leaf.u_shaded
  
  # 叶片尺度的净辐射更新方式
  # 参考Chen 2012年的聚集指数论文
  # return Rn_Leaf
  Rn_o, Rn_u, Rn_g
end
