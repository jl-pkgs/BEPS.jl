function netRadiation_c(shortRad_global, CosZs,
  temp_o, temp_u, temp_g,
  lai_o, lai_u, lai_os, lai_us, lai::Leaf, clumping, temp_air, rh,
  albedo_snow_v, albedo_snow_n,
  percentArea_snow_o, percentArea_snow_u, percent_snow_g,
  albedo_v_o, albedo_n_o, albedo_v_u, albedo_n_u, albedo_v_g, albedo_n_g,
  # netRad_o::TypeRef, netRad_u::TypeRef, netRad_g::TypeRef,
  netRadLeaf::Leaf, netShortRadLeaf::Leaf, Rnl_Leaf::Leaf, Ra::Radiation)

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


# @timeit_all 
function netRadiation_jl(Rs_global::FT, CosZs::FT,
  temp_o::FT, temp_u::FT, temp_g::FT,
  lai_o::FT, lai_u::FT, lai_os::FT, lai_us::FT,
  lai::Leaf,
  clumping::FT, temp_air::FT, rh::FT,
  albedo_snow_v::FT, albedo_snow_n::FT,
  percArea_snow_o::FT, percArea_snow_u::FT, perc_snow_g::FT,
  albedo_v_o::FT, albedo_n_o::FT, albedo_v_u::FT, albedo_n_u::FT,
  albedo_v_g::FT, albedo_n_g::FT,
  # Rn_o::FT, Rn_u::FT, Rn_g::FT,
  Rn_Leaf::Leaf,
  Rns_Leaf::Leaf,
  Rnl_Leaf::Leaf, Ra::Radiation)
  # Rnl_Leaf::Leaf = Leaf()

  # calculate albedo of canopy in this step
  albedo_v_os::FT = albedo_v_o * (1.0 - percArea_snow_o) + albedo_snow_v * percArea_snow_o  # visible, overstory
  albedo_n_os::FT = albedo_n_o * (1.0 - percArea_snow_o) + albedo_snow_n * percArea_snow_o  # near infrared
  albedo_v_us::FT = albedo_v_u * (1.0 - percArea_snow_u) + albedo_snow_v * percArea_snow_u  # understory
  albedo_n_us::FT = albedo_n_u * (1.0 - percArea_snow_u) + albedo_snow_n * percArea_snow_u

  albedo_o::FT = 0.5 * (albedo_v_os + albedo_n_os)
  albedo_u::FT = 0.5 * (albedo_v_us + albedo_n_us)

  # calculate albedo of ground in this step
  albedo_v_gs::FT = albedo_v_g * (1.0 - perc_snow_g) + albedo_snow_v * perc_snow_g
  albedo_n_gs::FT = albedo_n_g * (1.0 - perc_snow_g) + albedo_snow_n * perc_snow_g
  albedo_g::FT = 0.5 * (albedo_v_gs + albedo_n_gs)

  # separate global solar radiation into direct and diffuse one
  if (CosZs < 0.001)  # solar zenith angle small, all diffuse radiation
    ratio_cloud = 0.0
  else
    ratio_cloud = Rs_global / (1367 * CosZs)  # Luo2018, A4
  end

  if (ratio_cloud > 0.8)
    Ra.Rs_df = 0.13 * Rs_global  # Luo2018, A2
  else
    Ra.Rs_df = (0.943 + 0.734 * ratio_cloud - 4.9 * pow((ratio_cloud), 2) +
             1.796 * pow((ratio_cloud), 3) + 2.058 * pow((ratio_cloud), 4)) * Rs_global  # Luo2018, A2
  end

  Ra.Rs_df = clamp(Ra.Rs_df, 0.0, Rs_global)
  Ra.Rs_dir = Rs_global - Ra.Rs_df  # Luo2018, A3

  # fraction at each layer of canopy, direct and diffuse. use Leaf only lai here
  gap_o_dir::FT = exp(-0.5 * clumping * lai_o / CosZs)
  gap_u_dir::FT = exp(-0.5 * clumping * lai_u / CosZs)

  # indicators to describe leaf distribution angles in canopy. slightly related with LAI
  cosQ_o::FT = 0.537 + 0.025 * lai_o  # Luo2018, A10, a representative zenith angle for diffuse radiation transmission
  cosQ_u::FT = 0.537 + 0.025 * lai_u

  gap_o_df::FT = exp(-0.5 * clumping * lai_o / cosQ_o)
  gap_os_df::FT = exp(-0.5 * clumping * lai_os / cosQ_o)  # considering stem

  gap_u_df::FT = exp(-0.5 * clumping * lai_u / cosQ_u)
  gap_us_df::FT = exp(-0.5 * clumping * lai_us / cosQ_u)

  # emissivity of each part
  ea = cal_ea(temp_air, rh)
  emiss_air = 1.0 - exp(-(pow(ea * 10.0, (temp_air + 273.15) / 1200.0)))
  emiss_air = clamp(emiss_air, 0.7, 1.0)

  emiss_o = 0.98
  emiss_u = 0.98
  emiss_g = 0.96

  if Rs_global > zero && CosZs > zero
    # net short direct radiation on canopy and ground
    Ra.Rns_o_dir = Ra.Rs_dir * ((1.0 - albedo_o) - (1.0 - albedo_u) * gap_o_dir)  # dir into dif_under
    Ra.Rns_u_dir = Ra.Rs_dir * gap_o_dir * ((1.0 - albedo_u) - (1.0 - albedo_g) * gap_u_dir)
    Ra.Rns_g_dir = Ra.Rs_dir * gap_o_dir * gap_u_dir * (1.0 - albedo_g)
  else
    Ra.Rns_o_dir = 0.0
    Ra.Rns_u_dir = 0.0
    Ra.Rns_g_dir = 0.0
  end

  if Rs_global > zero && CosZs > zero
    # net short diffuse radiation on canopy and ground
    Ra.Rns_o_df = Ra.Rs_df * ((1.0 - albedo_o) - (1.0 - albedo_u) * gap_o_df) +
               0.21 * clumping * Ra.Rs_dir * (1.1 - 0.1 * lai_o) * exp(-CosZs)  # A8
    Ra.Rns_u_df = Ra.Rs_df * gap_o_df * ((1.0 - albedo_u) - (1.0 - albedo_g) * gap_u_df) +
               0.21 * clumping * Ra.Rs_dir * gap_o_dir * (1.1 - 0.1 * lai_u) * exp(-CosZs)  # A9
    Ra.Rns_g_df = Ra.Rs_df * gap_o_df * gap_u_df * (1.0 - albedo_g)
  else
    Ra.Rns_o_df = 0.0
    Ra.Rns_u_df = 0.0
    Ra.Rns_g_df = 0.0
  end

  # total net shortwave radiation at canopy level
  Rns_o = Ra.Rns_o_dir + Ra.Rns_o_df
  Rns_u = Ra.Rns_u_dir + Ra.Rns_u_df
  Rns_g = Ra.Rns_g_dir + Ra.Rns_g_df

  # 计算植被和地面的净长波辐射
  Rl_air = cal_Rln(emiss_air, temp_air)
  Rl_o = cal_Rln(emiss_o, temp_o)
  Rl_u = cal_Rln(emiss_u, temp_u)
  Rl_g = cal_Rln(emiss_g, temp_g)

  Rnl_o = (emiss_o * (Rl_air + Rl_u * (1.0 - gap_u_df) + Rl_g * gap_u_df) - 2 * Rl_o) *
          (1.0 - gap_o_df) +
          emiss_o * (1.0 - emiss_u) * (1.0 - gap_u_df) * (Rl_air * gap_o_df + Rl_o * (1.0 - gap_o_df))

  Rnl_u = (emiss_u * (Rl_air * gap_o_df + Rl_o * (1.0 - gap_o_df) + Rl_g) - 2 * Rl_u) * (1.0 - gap_u_df) +
          (1.0 - emiss_g) * ((Rl_air * gap_o_df + Rl_o * (1.0 - gap_o_df)) * gap_u_df + Rl_u * (1.0 - gap_u_df)) +
          emiss_u * (1.0 - emiss_o) * (Rl_u * (1.0 - gap_u_df) + Rl_g * gap_u_df) * (1.0 - gap_o_df)

  Rnl_g = emiss_g * ((Rl_air * gap_o_df + Rl_o * (1.0 - gap_o_df)) * gap_u_df + Rl_u * (1.0 - gap_u_df)) -
          Rl_g + (1.0 - emiss_u) * Rl_g * (1.0 - gap_u_df)

  # 计算植被和地面的总净辐射
  Rn_o = Rns_o + Rnl_o
  Rn_u = Rns_u + Rnl_u
  Rn_g = Rns_g + Rnl_g

  if Rs_global > zero && CosZs > zero # only happens in day time, when sun is out
    Rs_o_dir = 0.5 * Ra.Rs_dir / CosZs
    Rs_o_dir = min(Rs_o_dir, 0.7 * 1362)
    Rs_u_dir = Rs_o_dir

    Rs_o_df = (Ra.Rs_df - Ra.Rs_df * gap_os_df) / lai_os + 0.07 * Ra.Rs_dir * (1.1 - 0.1 * lai_os) * exp(-CosZs)
    Rs_u_df = (Ra.Rs_df * gap_o_df - Ra.Rs_df * gap_o_df * gap_us_df) / lai_us +
              0.05 * Ra.Rs_dir * gap_o_dir * (1.1 - 0.1 * lai_us) * exp(-CosZs)
  else
    Rs_o_dir = 0.0
    Rs_u_dir = 0.0
    Rs_o_df = 0.0
    Rs_u_df = 0.0
  end

  
  # overstorey sunlit leaves
  Rns_Leaf.o_sunlit = (Rs_o_dir + Rs_o_df) * (1.0 - albedo_o)
  Rnl_Leaf.o_sunlit = lai.o_sunlit > 0.0 ? Rnl_o / lai_os : Rnl_o

  # overstorey shaded leaf
  Rns_Leaf.o_shaded = Rs_o_df * (1.0 - albedo_o) # diffuse
  Rnl_Leaf.o_shaded = lai.o_shaded > 0.0 ? Rnl_o / lai_os : Rnl_o

  # understorey sunlit leaf
  Rns_Leaf.u_sunlit = (Rs_u_dir + Rs_u_df) * (1.0 - albedo_u)
  Rnl_Leaf.u_sunlit = lai.u_sunlit > 0.0 ? Rnl_u / lai_us : Rnl_u

  Rns_Leaf.u_shaded = Rs_u_df * (1.0 - albedo_u)
  Rnl_Leaf.u_shaded = lai.u_shaded > 0.0 ? Rnl_u / lai_us : Rnl_u

  Rn_Leaf.o_sunlit = Rns_Leaf.o_sunlit + Rnl_Leaf.o_sunlit
  Rn_Leaf.o_shaded = Rns_Leaf.o_shaded + Rnl_Leaf.o_shaded
  Rn_Leaf.u_sunlit = Rns_Leaf.u_sunlit + Rnl_Leaf.u_sunlit
  Rn_Leaf.u_shaded = Rns_Leaf.u_shaded + Rnl_Leaf.u_shaded

  # 叶片尺度的净辐射更新方式
  # 参考Chen 2012年的聚集指数论文
  Rn_o, Rn_u, Rn_g
end
