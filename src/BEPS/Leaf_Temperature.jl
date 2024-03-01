function Leaf_Temperature_jl(Tair::Float64, Δ::Float64, γ::Float64, VPD::Float64, cp::Float64,
  Gw::Float64, Gww::Float64, Gh::Float64, Xc_sl::Float64, radiation::Float64, constrain::Bool=true)

  p_star = (Gw + Gww * Xc_sl) / γ
  Tc = Tair + (radiation - VPD * ρₐ * cp * p_star) / (ρₐ * cp * (Gh + Δ * p_star))

  constrain && (Tc = clamp(Tc, Tair - 3.0, Tair + 5.0))
  return Tc
end


function Leaf_Temperatures_jl(Tair::Float64, slope::Float64, γ::Float64, 
  VPD_air::Float64, Cp_ca::Float64,
  Gw::Leaf, Gww::Leaf, Gh::Leaf,
  Xcs_o::Float64, Xcl_o::Float64, 
  Xcs_u::Float64, Xcl_u::Float64,
  radiation::Leaf, Tc::Leaf)

  args = (Tair, slope, γ, VPD_air, Cp_ca)
  Tc.o_sunlit = Leaf_Temperature_jl(args...,
    Gw.o_sunlit, Gww.o_sunlit, Gh.o_sunlit, Xcs_o + Xcl_o, radiation.o_sunlit)

  Tc.o_shaded = Leaf_Temperature_jl(args...,
    Gw.o_shaded, Gww.o_shaded, Gh.o_shaded, Xcs_o + Xcl_o, radiation.o_shaded)

  Tc.u_sunlit = Leaf_Temperature_jl(args...,
    Gw.u_sunlit, Gww.u_sunlit, Gh.u_sunlit, Xcs_u + Xcl_u, radiation.u_sunlit)

  Tc.u_shaded = Leaf_Temperature_jl(args...,
    Gw.u_shaded, Gww.u_shaded, Gh.u_shaded, Xcs_u + Xcl_u, radiation.u_shaded)
end
