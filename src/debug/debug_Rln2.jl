# cal_Rln_up = cal_Rln_down
function cal_Rln_net(L_up, L_down, Tair, ϵ, τ)
  Rl = cal_Rln(ϵ, Tair)
  ((L_up + L_down) * ϵ - 2 * Rl) * (1 - τ)
end

function cal_Rln_down(Rdown_i2, Rup_i1, Tair, ϵ, τd)
  Rl_i1 = cal_Rln(ϵ, Tair) # Rl
  ## 如果考虑反射，不考虑透射
  ρ = 1 - ϵ
  τl = 0
  # τd: 透射的部分
  # (1 - τd) * (1 - ϵ): 未透射的部分发射率
  Rdown_i2 * (τd + (1 - τd) * τl) + 
    Rup_i1 * (1 - τd) * ρ + 
    Rl_i1 * (1 - τd) # Eq 14.122
end

function cal_Rln_up(Rdown_i2, Rup_i1, Tair, ϵ, τd)
  Rl_i1 = cal_Rln(ϵ, Tair) # Rl
  ## 如果考虑反射，不考虑透射
  ρ = 1 - ϵ
  τl = 1 - ϵ
  ## 之前版本
  # ρ = 0
  # τl = 1 - ϵ
  Rup_i1 * (τd + (1 - τd) * τl) +
  Rdown_i2 * (1 - τd) * ρ +
  Rl_i1 * (1 - τd) # Eq 14.122
end


function debug_Rln(
  To::FT=25.0, Tu::FT=30.0, Tg::FT=35.0,
  lai_o::FT = 2.0, lai_u::FT = 1.0, clumping::FT=0.5, 
  Tair::FT=20.0, rh::FT=80.0)
  # Rnl_Leaf::Leaf = Leaf()
  # indicators to describe leaf distribution angles in canopy. slightly related with LAI
  cosQ_o::FT = 0.537 + 0.025 * lai_o  # Luo2018, A10, a representative zenith angle for diffuse radiation transmission
  cosQ_u::FT = 0.537 + 0.025 * lai_u

  τₒ::FT = exp(-0.5 * clumping * lai_o / cosQ_o)
  τᵤ::FT = exp(-0.5 * clumping * lai_u / cosQ_u)

  ea = cal_ea(Tair, rh)
  ϵ_air = 1.0 - exp(-(pow(ea * 10.0, (Tair + 273.15) / 1200.0)))
  ϵ_air = clamp(ϵ_air, 0.7, 1.0)

  ϵ_o = 0.98
  ϵ_u = 0.98
  ϵ_g = 0.96

  # 计算植被和地面的净长波辐射
  L_a = cal_Rln(ϵ_air, Tair)
  L_o = cal_Rln(ϵ_o, To)
  L_u = cal_Rln(ϵ_u, Tu)
  L_g = cal_Rln(ϵ_g, Tg)
  
  Lg_up = L_g
  Lu_up = cal_Rln_up(Lg_up, Lo_down, Tu, ϵ_u, τᵤ)
  Lo_up = cal_Rln_up(Lu_up, L_a, To, ϵ_o, τₒ)

  Lo_down = cal_Rln_down(L_a, Lu_up, To, ϵ_o, τₒ)
  Lu_down = cal_Rln_down(Lo_down, Lg_up, Tu, ϵ_u, τᵤ)
  
  Lg_down = cal_Rln_down(Lu_down, 0, Tg, ϵ_g, τᵤ)
  # Lu_up = cal_Rln_out(Lg_up, Tu, ϵ_u, τᵤ)
  # Lu_up = Lg_up * (1 - ϵ_u * (1 - tau_u)) + L_u * (1 - tau_u)
  # Lo_up = Lu_up * (1 - ϵ_o * (1 - τₒ)) + L_o * (1 - τₒ)
  # Lo_net = (L_a + Lu_up) * ϵ_o * (1 - τₒ) - 2 * L_o * (1 - τₒ)
  Lo_net = cal_Rln_net(L_a, Lu_up, To, ϵ_o, τₒ)
  Lu_net = cal_Rln_net(Lu_down, Lg_up, Tu, ϵ_u, τᵤ)
  Lg_net = cal_Rln_net(Lg_down, 0, Tu, ϵ_u, τᵤ)
  Lg_net = ϵ_g * Lu_down - L_g
  @show Lo_net, Lu_net, Lg_net
  @show Lo_net + Lu_net + Lg_net

  ## 另外一种算法
  Rnl_o = (ϵ_o * (L_a + L_u * (1.0 - τᵤ) + L_g * τᵤ) - 2 * L_o) * (1.0 - τₒ) +
          ϵ_o * (1.0 - ϵ_u) * (1.0 - τᵤ) * (L_a * τₒ + L_o * (1.0 - τₒ))

  Rnl_u = (ϵ_u * (L_a * τₒ + L_o * (1.0 - τₒ) + L_g) - 2 * L_u) * (1.0 - τᵤ) +
          (1.0 - ϵ_g) * ((L_a * τₒ + L_o * (1.0 - τₒ)) * τᵤ + L_u * (1.0 - τᵤ)) +
          ϵ_u * (1.0 - ϵ_o) * (L_u * (1.0 - τᵤ) + L_g * τᵤ) * (1.0 - τₒ)

  # 这里它考了向后的反射
  Rnl_g = ϵ_g * ((L_a * τₒ + L_o * (1.0 - τₒ)) * τᵤ + L_u * (1.0 - τᵤ)) -
          L_g + L_g * (1.0 - ϵ_u) * (1.0 - τᵤ)

  @show (; Rnl_o, Rnl_u, Rnl_g)
  Rnl_o + Rnl_u + Rnl_g
end

export debug_Rln
