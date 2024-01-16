```julia
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
```
