"""
## Arguments
- `T`                 : temperature of o, u, g
- `α_v, α_n`          : albedo of visible, near infrared
- `percentArea_snow_o`: percentage of snow on overstorey (by area)
- `percentArea_snow_u`: percentage of snow on understorey (by area)
- `percent_snow_g`    : percentage of snow on ground (by mass)
"""
function netRadiation_jl(Rs_global::FT, CosZs::FT,
  T::Layer3{FT}, 
  lai_o::FT, lai_u::FT, lai_os::FT, lai_us::FT,
  lai::Leaf,
  Ω::FT, Tair::FT, RH::FT, LR::FT,
  α_snow_v::FT, α_snow_n::FT, α_v::Layer3{FT}, α_n::Layer3{FT},
  percArea_snow_o::FT, percArea_snow_u::FT, perc_snow_g::FT,
  Rn_Leaf::Leaf, Rns_Leaf::Leaf, Rnl_Leaf::Leaf, Ra::Radiation)
  
  # calculate α of canopy in this step
  α_v_os::FT = α_v.o * (1.0 - percArea_snow_o) + α_snow_v * percArea_snow_o  # visible, overstory
  α_n_os::FT = α_n.o * (1.0 - percArea_snow_o) + α_snow_n * percArea_snow_o  # near infrared
  α_v_us::FT = α_v.u * (1.0 - percArea_snow_u) + α_snow_v * percArea_snow_u  # understory
  α_n_us::FT = α_n.u * (1.0 - percArea_snow_u) + α_snow_n * percArea_snow_u

  α_o::FT = 0.5 * (α_v_os + α_n_os)
  α_u::FT = 0.5 * (α_v_us + α_n_us)

  # calculate α of ground in this step
  α_v_gs::FT = α_v.g * (1.0 - perc_snow_g) + α_snow_v * perc_snow_g
  α_n_gs::FT = α_n.g * (1.0 - perc_snow_g) + α_snow_n * perc_snow_g
  α_g::FT = 0.5 * (α_v_gs + α_n_gs)

  # separate global solar radiation into direct and diffuse one
  # solar zenith angle small, all diffuse radiation
  ratio_cloud = (CosZs < 0.001) ? 0.0 : Rs_global / (1367 * CosZs)  # Luo2018, A4
  
  if (ratio_cloud > 0.8)
    Ra.Rs_df = 0.13 * Rs_global  # Luo2018, A2
  else
    Ra.Rs_df = (0.943 + 0.734 * ratio_cloud - 4.9 * pow((ratio_cloud), 2) +
             1.796 * pow((ratio_cloud), 3) + 2.058 * pow((ratio_cloud), 4)) * Rs_global  # Luo2018, A2
  end

  Ra.Rs_df = clamp(Ra.Rs_df, 0.0, Rs_global)
  Ra.Rs_dir = Rs_global - Ra.Rs_df  # Luo2018, A3

  # fraction at each layer of canopy, direct and diffuse. use Leaf only lai here
  τ_o_dir::FT = exp(-0.5 * Ω * lai_o / CosZs)
  τ_u_dir::FT = exp(-0.5 * Ω * lai_u / CosZs)

  # indicators to describe leaf distribution angles in canopy. slightly related with LAI
  cosQ_o::FT = 0.537 + 0.025 * lai_o  # Luo2018, A10, a representative zenith angle for diffuse radiation transmission
  cosQ_u::FT = 0.537 + 0.025 * lai_u

  τ_o_df::FT = exp(-0.5 * Ω * lai_o / cosQ_o)
  τ_os_df::FT = exp(-0.5 * Ω * lai_os / cosQ_o)  # considering stem

  τ_u_df::FT = exp(-0.5 * Ω * lai_u / cosQ_u)
  τ_us_df::FT = exp(-0.5 * Ω * lai_us / cosQ_u)

  # net short direct radiation on canopy and ground
  if Rs_global > 0.0 && CosZs > 0.0
    Ra.Rns_o_dir = Ra.Rs_dir * ((1.0 - α_o) - (1.0 - α_u) * τ_o_dir)  # dir into dif_under
    Ra.Rns_u_dir = Ra.Rs_dir * τ_o_dir * ((1.0 - α_u) - (1.0 - α_g) * τ_u_dir)
    Ra.Rns_g_dir = Ra.Rs_dir * τ_o_dir * τ_u_dir * (1.0 - α_g)
  else
    Ra.Rns_o_dir = 0.0
    Ra.Rns_u_dir = 0.0
    Ra.Rns_g_dir = 0.0
  end

  # net short diffuse radiation on canopy and ground
  if Rs_global > 0.0 && CosZs > 0.0
    Ra.Rns_o_df = Ra.Rs_df * ((1.0 - α_o) - (1.0 - α_u) * τ_o_df) +
               0.21 * Ω * Ra.Rs_dir * (1.1 - 0.1 * lai_o) * exp(-CosZs)  # A8
    Ra.Rns_u_df = Ra.Rs_df * τ_o_df * ((1.0 - α_u) - (1.0 - α_g) * τ_u_df) +
               0.21 * Ω * Ra.Rs_dir * τ_o_dir * (1.1 - 0.1 * lai_u) * exp(-CosZs)  # A9
    Ra.Rns_g_df = Ra.Rs_df * τ_o_df * τ_u_df * (1.0 - α_g)
  else
    Ra.Rns_o_df = 0.0
    Ra.Rns_u_df = 0.0
    Ra.Rns_g_df = 0.0
  end

  # total net shortwave radiation at canopy level
  Rns_o = Ra.Rns_o_dir + Ra.Rns_o_df
  Rns_u = Ra.Rns_u_dir + Ra.Rns_u_df
  Rns_g = Ra.Rns_g_dir + Ra.Rns_g_df

  Rnl_o, Rnl_u, Rnl_g = cal_Rln_Longwave(Tair, RH, T, lai_o, lai_u, Ω, LR)

  # 计算植被和地面的总净辐射
  Rn_o = Rns_o + Rnl_o
  Rn_u = Rns_u + Rnl_u
  Rn_g = Rns_g + Rnl_g

  if Rs_global > 0.0 && CosZs > 0.0 # only happens in day time, when sun is out
    Rs_o_dir = 0.5 * Ra.Rs_dir / CosZs
    Rs_o_dir = min(Rs_o_dir, 0.7 * 1362)
    Rs_u_dir = Rs_o_dir

    Rs_o_df = (Ra.Rs_df - Ra.Rs_df * τ_os_df) / lai_os + 0.07 * Ra.Rs_dir * (1.1 - 0.1 * lai_os) * exp(-CosZs)
    Rs_u_df = (Ra.Rs_df * τ_o_df - Ra.Rs_df * τ_o_df * τ_us_df) / lai_us +
              0.05 * Ra.Rs_dir * τ_o_dir * (1.1 - 0.1 * lai_us) * exp(-CosZs)
  else
    Rs_o_dir = 0.0
    Rs_u_dir = 0.0
    Rs_o_df = 0.0
    Rs_u_df = 0.0
  end

  # overstorey sunlit leaves
  Rns_Leaf.o_sunlit = (Rs_o_dir + Rs_o_df) * (1.0 - α_o)
  Rnl_Leaf.o_sunlit = lai.o_sunlit > 0.0 ? Rnl_o / lai_os : Rnl_o

  # overstorey shaded leaf
  Rns_Leaf.o_shaded = Rs_o_df * (1.0 - α_o) # diffuse
  Rnl_Leaf.o_shaded = lai.o_shaded > 0.0 ? Rnl_o / lai_os : Rnl_o

  # understorey sunlit leaf
  Rns_Leaf.u_sunlit = (Rs_u_dir + Rs_u_df) * (1.0 - α_u)
  Rnl_Leaf.u_sunlit = lai.u_sunlit > 0.0 ? Rnl_u / lai_us : Rnl_u

  Rns_Leaf.u_shaded = Rs_u_df * (1.0 - α_u)
  Rnl_Leaf.u_shaded = lai.u_shaded > 0.0 ? Rnl_u / lai_us : Rnl_u

  Rn_Leaf.o_sunlit = Rns_Leaf.o_sunlit + Rnl_Leaf.o_sunlit
  Rn_Leaf.o_shaded = Rns_Leaf.o_shaded + Rnl_Leaf.o_shaded
  Rn_Leaf.u_sunlit = Rns_Leaf.u_sunlit + Rnl_Leaf.u_sunlit
  Rn_Leaf.u_shaded = Rns_Leaf.u_shaded + Rnl_Leaf.u_shaded

  # 叶片尺度的净辐射更新方式
  # 参考Chen 2012年的聚集指数论文
  Rn_o, Rn_u, Rn_g
end


function cal_Rln_Longwave(Tair::FT, RH::FT, T::Layer3{FT}, 
  lai_o::FT, lai_u::FT, Ω::FT, LR::FT) where {FT<:Real}

  # indicators to describe leaf distribution angles in canopy. slightly related with LAI
  cosQ_o::FT = 0.537 + 0.025 * lai_o  # Luo2018, A10, a representative zenith angle for diffuse radiation transmission
  cosQ_u::FT = 0.537 + 0.025 * lai_u

  τ_o_df::FT = exp(-0.5 * Ω * lai_o / cosQ_o)
  τ_u_df::FT = exp(-0.5 * Ω * lai_u / cosQ_u)

  # ϵ of each part
  ea = cal_ea(Tair, RH)
  ϵ_air = 1.0 - exp(-(pow(ea * 10.0, (Tair + 273.15) / 1200.0)))
  ϵ_air = clamp(ϵ_air, 0.7, 1.0)

  ϵ_o = 0.98
  ϵ_u = 0.98
  ϵ_g = 0.96

  # 计算植被和地面的净长波辐射
  Rl_air = isfinite(LR) ? LR : cal_Rln(ϵ_air, Tair)
  Rl_o = cal_Rln(ϵ_o, T.o)
  Rl_u = cal_Rln(ϵ_u, T.u)
  Rl_g = cal_Rln(ϵ_g, T.g)

  Rnl_o = (ϵ_o * (Rl_air + Rl_u * (1.0 - τ_u_df) + Rl_g * τ_u_df) - 2 * Rl_o) *
          (1.0 - τ_o_df) +
          ϵ_o * (1.0 - ϵ_u) * (1.0 - τ_u_df) * (Rl_air * τ_o_df + Rl_o * (1.0 - τ_o_df))

  Rnl_u = (ϵ_u * (Rl_air * τ_o_df + Rl_o * (1.0 - τ_o_df) + Rl_g) - 2 * Rl_u) * (1.0 - τ_u_df) +
          (1.0 - ϵ_g) * ((Rl_air * τ_o_df + Rl_o * (1.0 - τ_o_df)) * τ_u_df + Rl_u * (1.0 - τ_u_df)) +
          ϵ_u * (1.0 - ϵ_o) * (Rl_u * (1.0 - τ_u_df) + Rl_g * τ_u_df) * (1.0 - τ_o_df)

  Rnl_g = ϵ_g * ((Rl_air * τ_o_df + Rl_o * (1.0 - τ_o_df)) * τ_u_df + Rl_u * (1.0 - τ_u_df)) -
          Rl_g + (1.0 - ϵ_u) * Rl_g * (1.0 - τ_u_df)
  
  Rnl_o, Rnl_u, Rnl_g
end
