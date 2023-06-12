pow = ^

# kPa deg-1
cal_slope(Ta::Real) = 2503.0 / pow((Ta + 237.3), 2) * exp(17.27 * Ta / (Ta + 237.3))


# kPa
cal_es(Ta::Real) = 0.61078 * exp(17.3 * Ta / (237.3 + Ta))

cal_ea(Ta::Real, RH::Real) = cal_es(Ta) * RH / 100

cal_lambda(Ta::Real) = (2.501 - 0.00237 * Ta) * 1000000

ea2q(ea::Real, Pa::Real=101.35) = 0.622 * ea / (Pa - 0.378 * ea)

function RH2q(Ta::Real, RH::Real)
  es = cal_es(Ta)
  ea = es * RH / 100
  ea2q(ea)
end


cal_cp(q::Real) = 1004.65 * (1 + 0.84 * q)

function cal_cp(Ta::Real, RH::Real)
  q = RH2q(Ta, RH)
  cal_cp(q)
end 


function meteo_pack_jl(Ta::Real, RH::Real)
  rho_a = 1.292 # ρ_air, kg/m3

  es = cal_es(Ta)
  ea = es * RH / 100
  vpd = es - ea

  q = ea2q(ea)
  cp = 1004.65 * (1 + 0.84 * q)
  
  slope = cal_slope(Ta) # slope of es
  psy = 0.066            # psychrometer γ, kPa/C, 
  # lambda = cal_lambda(Ta) # J kg-1
  # psy = cp * 101.13 / (0.622 * lambda)
  (; rho_a, cp, vpd, slope, psy, es, ea, q)
end


function update_Gw!(Gw::Leaf, Gs_new::Leaf, Ga_o, Ga_u, Gb_o, Gb_u)
  Gw.o_sunlit = 1.0 / (1.0 / Ga_o + 1.0 / Gb_o + 1.0 / Gs_new.o_sunlit)  #conductance for water
  Gw.o_shaded = 1.0 / (1.0 / Ga_o + 1.0 / Gb_o + 1.0 / Gs_new.o_shaded)
  Gw.u_sunlit = 1.0 / (1.0 / Ga_u + 1.0 / Gb_u + 1.0 / Gs_new.u_sunlit)
  Gw.u_shaded = 1.0 / (1.0 / Ga_u + 1.0 / Gb_u + 1.0 / Gs_new.u_shaded)
end

function update_Gc!(Gc::Leaf, Gs_new::Leaf, Ga_o, Ga_u, Gb_o, Gb_u)
  Gc.o_sunlit = 1.0 / (1.0 / Ga_o + 1.4 / Gb_o + 1.6 / Gs_new.o_sunlit) # conductance for CO2
  Gc.o_shaded = 1.0 / (1.0 / Ga_o + 1.4 / Gb_o + 1.6 / Gs_new.o_shaded)
  Gc.u_sunlit = 1.0 / (1.0 / Ga_u + 1.4 / Gb_u + 1.6 / Gs_new.u_sunlit)
  Gc.u_shaded = 1.0 / (1.0 / Ga_u + 1.4 / Gb_u + 1.6 / Gs_new.u_shaded)
end


function update_Gw!(Gw::Leaf, Gs_new::LeafRef, Ga_o, Ga_u, Gb_o, Gb_u)
  Gw.o_sunlit = 1.0 / (1.0 / Ga_o + 1.0 / Gb_o + 1.0 / Gs_new.o_sunlit[])  #conductance for water
  Gw.o_shaded = 1.0 / (1.0 / Ga_o + 1.0 / Gb_o + 1.0 / Gs_new.o_shaded[])
  Gw.u_sunlit = 1.0 / (1.0 / Ga_u + 1.0 / Gb_u + 1.0 / Gs_new.u_sunlit[])
  Gw.u_shaded = 1.0 / (1.0 / Ga_u + 1.0 / Gb_u + 1.0 / Gs_new.u_shaded[])
end

function update_Gc!(Gc::Leaf, Gs_new::LeafRef, Ga_o, Ga_u, Gb_o, Gb_u)
  Gc.o_sunlit = 1.0 / (1.0 / Ga_o + 1.4 / Gb_o + 1.6 / Gs_new.o_sunlit[]) # conductance for CO2
  Gc.o_shaded = 1.0 / (1.0 / Ga_o + 1.4 / Gb_o + 1.6 / Gs_new.o_shaded[])
  Gc.u_sunlit = 1.0 / (1.0 / Ga_u + 1.4 / Gb_u + 1.6 / Gs_new.u_sunlit[])
  Gc.u_shaded = 1.0 / (1.0 / Ga_u + 1.4 / Gb_u + 1.6 / Gs_new.u_shaded[])
end


export meteo_pack_jl, cal_slope, cal_es
