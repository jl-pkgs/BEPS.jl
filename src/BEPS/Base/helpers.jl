# pow = ^
pow(x::FT, y)::FT = x^y

function blackbody(T::FT)
  sb_constant = 5.67 / 100000000    # stephen-boltzman constant
  sb_constant * (T + 273.15)^4
end

function cal_Rln(emiss::FT, T::FT)
  emiss * blackbody(T)
end

# kPa deg-1
cal_slope(Ta::FT)::FT = 2503.0 / pow((Ta + 237.3), 2) * exp(17.27 * Ta / (Ta + 237.3))

# kPa
cal_es(Ta::FT)::FT = 0.61078 * exp(17.3 * Ta / (237.3 + Ta))

cal_ea(Ta::FT, RH::FT)::FT = cal_es(Ta) * RH / 100

cal_lambda(Ta::FT)::FT = (2.501 - 0.00237 * Ta) * 1000000

ea2q(ea::FT, Pa::FT=101.35)::FT = 0.622 * ea / (Pa - 0.378 * ea)

function RH2q(Ta::FT, RH::FT)::FT
  es = cal_es(Ta)
  ea = es * RH / 100
  ea2q(ea)
end


cal_cp(q::FT)::FT = 1004.65 * (1 + 0.84 * q)

function cal_cp(Ta::FT, RH::FT)::FT
  q = RH2q(Ta, RH)
  cal_cp(q)
end 


function meteo_pack_jl(Ta::FT, RH::FT)
  rho_a = 1.292 # ρ_air, kg/m3

  es = cal_es(Ta)
  ea = es * RH / 100
  VPD = es - ea

  q = ea2q(ea)
  cp = cal_cp(q)
  
  slope = cal_slope(Ta) # slope of es
  gamma = 0.066            # γ γ, kPa/C, 
  # lambda = cal_lambda(Ta) # J kg-1
  # psy = cp * 101.13 / (0.622 * lambda)
  (; rho_a, cp, VPD, slope, gamma, es, ea, q)
end



export pow, cal_ea, cal_es, 
  cal_slope, cal_lambda, cal_cp, ea2q, RH2q,
  cal_Rln, 
  meteo_pack_jl, blackbody
