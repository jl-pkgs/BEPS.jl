# L_down = L_down_2 * (τ + (1 - τ) * τ_l) + L_up * (1 - τ) * ρ + cal_Rln(ϵ, Tair) * (1 - τ)
# L_up_2 = L_up * (τ + (1 - τ) * t_l) + L_down_2 * (1 - τ) * ρ + cal_Rln(ϵ, Tair) * (1 - τ)
# # 这一层的净辐射为
# Rnl = L_down_2 + L_up - L_down - L_up_2

# cal_Rln_up = cal_Rln_down
function cal_Rln_net(L_up, L_down, Tair, ϵ, τ)
  Rl = cal_Rln(ϵ, Tair)
  ((L_up + L_down) * ϵ - 2 * Rl) * (1 - τ)
end

function Rl_down()
  # τ[i+1] * (1 - ϵ[i]) + ϵ[i] * (1 - τ[i]``
  # 全部都是i+1
  # I[i+1] * (1 - ϵ[i+1] * (1 - τ[i+1])) + cal_Rln(ϵ[i+1], T[i+1]) * (1 - τ[i+1])
  # I[i+1] * (1 - ϵ[i+1] * (1 - τ[i+1])) + cal_Rln(ϵ[i+1], T[i+1]) * (1 - τ[i+1])
  L_sky * (1 - ϵ_o * (1 - τ_o)) + L_o * (1 - τ_o)
end
# blackbody(epsilon, T) = sigma * epsilon * T^4;

function cal_Rln_out(Rl_in, Tair, ϵ, τ)
  Rl = cal_Rln(ϵ, Tair)
  Rl_in * ((1 - ϵ) * (1 - τ) + τ) + Rl * (1 - τ)
  # Rl_in * (1 - τ)
  # Rl_in * (1 - ϵ * (1 - τ)) + Rl * (1 - τ)
end

function debug_Rln(
  To::FT=25.0, Tu::FT=30.0, Tg::FT=35.0,
  lai_o::FT=2.0, lai_u::FT=1.0, clumping::FT=0.5,
  Tair::FT=20.0, RH::FT=80.0)
  # Rnl_Leaf::Leaf = Leaf()
  # indicators to describe leaf distribution angles in canopy. slightly related with LAI
  cosQ_o::FT = 0.537 + 0.025 * lai_o  # Luo2018, A10, a representative zenith angle for diffuse radiation transmission
  cosQ_u::FT = 0.537 + 0.025 * lai_u

  τ_o::FT = exp(-0.5 * clumping * lai_o / cosQ_o) # 
  τ_u::FT = exp(-0.5 * clumping * lai_u / cosQ_u)
  τ_g = 1
  
  ea = cal_ea(Tair, RH)
  ϵ_air = 1.0 - exp(-(pow(ea * 10.0, (Tair + 273.15) / 1200.0)))
  ϵ_air = clamp(ϵ_air, 0.7, 1.0)

  ϵ_o = 0.98
  ϵ_u = 0.98
  ϵ_g = 0.96

  L_sky = cal_Rln(ϵ_air, Tair)
  L_o = cal_Rln(ϵ_o, To)
  L_u = cal_Rln(ϵ_u, Tu)
  L_g = cal_Rln(ϵ_g, Tg)

  # 反射率
  ρ_o = 1 - ϵ_o
  ρ_u = 1 - ϵ_u
  ρ_g = 1 - ϵ_g

  Rl_up_u = L_u * (1.0 - τ_u) + L_g * τ_u         # u向上的部分
  # Rl_up_o = L_o * (1.0 - τ_o) + Rl_up_u * τ_o

  Rl_down_o = (L_sky * τ_o + L_o * (1.0 - τ_o))       # o向下的部分
  Rl_down_u = (Rl_down_o * τ_u + L_u * (1.0 - τ_u)) # u向下的部分

  Lf_u2o = Rl_down_o * (1.0 - τ_u) * ρ_u # Lf代表反射的部分
  Lf_o2u = Rl_up_u * (1.0 - τ_o) * ρ_o
  Lf_g2u = Rl_down_u * (1.0 - τ_g) * ρ_g
  Lf_u2g = L_g * (1.0 - τ_u) * ρ_u
  
  # Lf_u2g, Lf_g2u, Lf_u2o, Lf_o2u = 0, 0, 0, 0
  
  # (1.0 - τ_o)：能量的吸收率
  # (1.0 - ϵ_o): （发射率=吸收率），未吸收的部分认为是反射（吸收能量反射的部分）
  # (1.0 - τ_o) * (1.0 - ϵ_o): 吸收的能量一部分反射出去，只考虑向后的反射  
  Rnl_o = (ϵ_o * (L_sky + Rl_up_u + Lf_u2o) - 2 * L_o) * (1.0 - τ_o)
  Rnl_u = (ϵ_u * (Rl_down_o + L_g + Lf_o2u + Lf_g2u) - 2 * L_u) * (1.0 - τ_u)
  Rnl_g = (ϵ_g * (Rl_down_u + Lf_u2g) - L_g) * (1.0 - τ_u)

  # display((; Rnl_o, Rnl_u, Rnl_g))
  # Rnl_o + Rnl_u + Rnl_g
  # -97.56872055403599
  Rnl_o, Rnl_u, Rnl_g
end

# # ϵ_o, ϵ_u, ϵ_g = ϵ_o, ϵ_u, ϵ_g
# Lo_down = cal_Rln_out(L_sky, To, ϵ_o, τ_o)
# # Lo_down = L_sky * (1 - ϵ_o * (1 - τ_o)) + L_o * (1 - τ_o)
# Rl_down_o = (L_sky * τ_o + L_o * (1.0 - τ_o))
# Lu_down = cal_Rln_out(Lo_down, Tu, ϵ_u, τ_u)
# # Lu_down = Lo_down * (1 - ϵ_u * (1 - tau_u)) + L_u * (1 - tau_u)
# Lg_up = L_g + (1 - ϵ_g) * Lu_down
# # Lg_up = L_g
# Lu_up = cal_Rln_out(Lg_up, Tu, ϵ_u, τ_u)
# # Lu_up = Lg_up * (1 - ϵ_u * (1 - tau_u)) + L_u * (1 - tau_u)
# # Lo_up = Lu_up * (1 - ϵ_o * (1 - τ_o)) + L_o * (1 - τ_o)
# # Lo_net = (L_sky + Lu_up) * ϵ_o * (1 - τ_o) - 2 * L_o * (1 - τ_o)

# Lo_net = cal_Rln_net(Lu_up, L_sky, To, ϵ_o, τ_o)
# Lu_net = cal_Rln_net(Lg_up, Lo_down, Tu, ϵ_u, τ_u)
# Lg_net = cal_Rln_net(Lg_up, Lo_down, Tu, ϵ_u, τ_u)
# Lg_net = ϵ_g * Lu_down - L_g
# @show Lo_net, Lu_net, Lg_net
# @show Lo_net + Lu_net + Lg_net

## 另外一种算法
# L_sky * τ_o: 到达i+2->i的长波辐射
# L_o * (1.0 - τ_o): 到达i+1->i的长波辐射

# @show τ_o, τ_u
# k = (1 - τ_o) * ρ_o * (1 - τ_u) * ρ_u


# x * τ_o
# x * (1 - τ_o) * ϵ_o
# x * (1 - τ_o) * ρ_o # 向下反射的部分

# x * (1 - τ_u) * ρ_u # 再次反射 ↑

# x * τ_o 
# k = (1 - τ_o) * ρ_o * (1 - τ_u) * ρ_u
# x * k # 再次反射 ↑
# x (k^0 + k^1 + k^2)



export debug_Rln
