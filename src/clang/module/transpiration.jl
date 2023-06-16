function transpiration_c(T_leaf::Leaf, Ta::Float64, RH::Float64, Gtrans::Leaf, lai::Leaf)

  trans_o = init_dbl()
  trans_u = init_dbl()

  ccall((:transpiration, libbeps), Cvoid,
    (Leaf, Cdouble, Cdouble, Leaf, Leaf, Ptr{Cdouble}, Ptr{Cdouble}),
    T_leaf, Ta, RH, Gtrans, lai, trans_o, trans_u)
  trans_o[], trans_u[]
end


function transpiration_jl(T_leaf::Leaf, Ta::Float64, RH::Float64, Gtrans::Leaf, lai::Leaf)
  T = Leaf() # transpiration
  
  met = meteo_pack_jl(Ta, RH)
  ρₐ = rho_a
  cp = met.cp     # specific heat of moist air above canopy
  VPD = met.VPD    # VPD in air
  Δ = met.slope
  γ = met.gamma  # psychrometer constant
  λ = cal_lambda(Ta)

  T.o_sunlit = (VPD + Δ * (T_leaf.o_sunlit - Ta)) * ρₐ * cp * Gtrans.o_sunlit / γ
  T.o_shaded = (VPD + Δ * (T_leaf.o_shaded - Ta)) * ρₐ * cp * Gtrans.o_shaded / γ
  T.u_sunlit = (VPD + Δ * (T_leaf.u_sunlit - Ta)) * ρₐ * cp * Gtrans.u_sunlit / γ
  T.u_shaded = (VPD + Δ * (T_leaf.u_shaded - Ta)) * ρₐ * cp * Gtrans.u_shaded / γ

  trans_o = 1.0 / λ * (T.o_sunlit * lai.o_sunlit + T.o_shaded * lai.o_shaded)
  trans_u = 1.0 / λ * (T.u_sunlit * lai.u_sunlit + T.u_shaded * lai.u_shaded)

  trans_o, trans_u
end
