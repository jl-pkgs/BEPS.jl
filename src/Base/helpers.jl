# pow = ^
pow(x::FT, y)::FT = x^y

function blackbody(T::FT)
  σ = 5.67 / 100000000    # stephen-boltzman constant
  σ * (T + 273.15)^4
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
  ρₐ::FT = 1.292 # ρₐir, kg/m3

  es::FT = cal_es(Ta)
  ea::FT = es * RH / 100
  VPD::FT = es - ea

  q::FT = ea2q(ea)
  cp::FT = cal_cp(q)

  λ::FT = cal_lambda(Ta)
  Δ::FT = cal_slope(Ta) # slope of es
  γ::FT = 0.066         # kPa/K, 
  # lambda = cal_lambda(Ta) # J kg-1
  # psy = cp * 101.13 / (0.622 * lambda)  
  (; ρₐ, cp, VPD, λ, Δ, γ, es, ea, q)
end



export pow, cal_ea, cal_es,
  cal_slope, cal_lambda, cal_cp, ea2q, RH2q,
  cal_Rln,
  meteo_pack_jl, blackbody
