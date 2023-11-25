function cal_Rln_out(Rl_in, Tair, ϵ, τ)
  Rl = cal_Rln(ϵ, Tair)
  Rl_in * (1 - ϵ * (1 - τ)) + Rl * (1 - τ)
end

# cal_Rln_up = cal_Rln_down
function cal_Rln_net(L_up, L_down, Tair, ϵ, τ)
  Rl = cal_Rln(ϵ, Tair)
  ((L_up + L_down) * ϵ - 2 * Rl) * (1 - τ)
end

function Rl_down()
  # τ[i+1] * (1 - ϵ[i]) + ϵ[i] * (1 - τ[i]
  # 全部都是i+1
  # I[i+1] * (1 - ϵ[i+1] * (1 - τ[i+1])) + cal_Rln(ϵ[i+1], T[i+1]) * (1 - τ[i+1])
  # I[i+1] * (1 - ϵ[i+1] * (1 - τ[i+1])) + cal_Rln(ϵ[i+1], T[i+1]) * (1 - τ[i+1])
  L_a * (1 - ϵ_o * (1 - τₒ)) + L_o * (1 - τₒ)
end
# blackbody(epsilon, T) = sigma * epsilon * T^4;


function debug_Rln(
  To::FT=25.0, Tu::FT=30.0, Tg::FT=35.0,
  lai_o::FT=2.0, lai_u::FT=1.0, clumping::FT=0.5,
  Tair::FT=20.0, rh::FT=80.0)
  # Rnl_Leaf::Leaf = Leaf()
  # indicators to describe leaf distribution angles in canopy. slightly related with LAI
  cosQ_o::FT = 0.537 + 0.025 * lai_o  # Luo2018, A10, a representative zenith angle for diffuse radiation transmission
  cosQ_u::FT = 0.537 + 0.025 * lai_u

  τₒ::FT = exp(-0.5 * clumping * lai_o / cosQ_o)
  τᵤ::FT = exp(-0.5 * clumping * lai_u / cosQ_u)

  # gap_os_df::FT = exp(-0.5 * clumping * lai_os / cosQ_o)  # considering stem
  # gap_us_df::FT = exp(-0.5 * clumping * lai_us / cosQ_u)
  # ϵivity of each part
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

  # ϵ_o, ϵ_u, ϵ_g = ϵ_o, ϵ_u, ϵ_g
  Lo_down = cal_Rln_out(L_a, To, ϵ_o, τₒ)
  # Lo_down = L_a * (1 - ϵ_o * (1 - τₒ)) + L_o * (1 - τₒ)

  Lu_down = cal_Rln_out(Lo_down, Tu, ϵ_u, τᵤ)
  # Lu_down = Lo_down * (1 - ϵ_u * (1 - tau_u)) + L_u * (1 - tau_u)
  Lg_up = L_g + (1 - ϵ_g) * Lu_down
  # Lg_up = L_g

  Lu_up = cal_Rln_out(Lg_up, Tu, ϵ_u, τᵤ)
  # Lu_up = Lg_up * (1 - ϵ_u * (1 - tau_u)) + L_u * (1 - tau_u)
  # Lo_up = Lu_up * (1 - ϵ_o * (1 - τₒ)) + L_o * (1 - τₒ)
  # Lo_net = (L_a + Lu_up) * ϵ_o * (1 - τₒ) - 2 * L_o * (1 - τₒ)
  Lo_net = cal_Rln_net(Lu_up, L_a, To, ϵ_o, τₒ)
  Lu_net = cal_Rln_net(Lg_up, Lo_down, Tu, ϵ_u, τᵤ)
  Lg_net = cal_Rln_net(Lg_up, Lo_down, Tu, ϵ_u, τᵤ)
  Lg_net = ϵ_g * Lu_down - L_g
  @show Lo_net, Lu_net, Lg_net
  @show Lo_net + Lu_net + Lg_net

  ## 另外一种算法
  # L_a * τₒ: 到达i+2->i的长波辐射
  # L_o * (1.0 - τₒ): 到达i+1->i的长波辐射
Rl_down_o = (L_a * τₒ + L_o * (1.0 - τₒ))       # o向下的部分
Rl_down_u = (Rl_down_o * τᵤ + L_u * (1.0 - τᵤ)) # u向下的部分
Rl_up_u   = L_u * (1.0 - τᵤ) + L_g * τᵤ           # u向上的部分

# 反射率
ρ_o = 1 - ϵ_o
ρ_u = 1 - ϵ_u
ρ_g = 1 - ϵ_g

# (1.0 - τₒ)：能量的吸收率
# (1.0 - ϵ_o): （发射率=吸收率），未吸收的部分认为是反射（吸收能量反射的部分）
# (1.0 - τₒ) * (1.0 - ϵ_o): 吸收的能量一部分反射出去，只考虑向后的反射
Rnl_o = (ϵ_o * (L_a + Rl_up_u) - 2 * L_o) * (1.0 - τₒ) +
        ϵ_o * Rl_down_o * (1.0 - τᵤ) * ρ_o  # understory反射的部分

Rnl_u = (ϵ_u * (Rl_down_o + L_g) - 2 * L_u) * (1.0 - τᵤ) +
        ρ_g * Rl_down_u +                 # 地表反射的部分, miss ϵ_u
        ϵ_u * Rl_up_u * (1.0 - τₒ) * ρ_o  # o反射的部分

Rnl_g = ϵ_g * Rl_down_u - L_g +
        L_g * (1.0 - τᵤ) * ρ_o # u反射的部分, miss a ϵ_g

  display((; Rnl_o, Rnl_u, Rnl_g))
  Rnl_o + Rnl_u + Rnl_g
  # -97.56872055403599
end

export debug_Rln
