# pow = ^
# pow(x::FT, y)::FT where {FT<:Real} = x^y
pow(x, y) = x^y

function blackbody(T::FT) where {FT<:Real}
  σ = 5.67 / 100000000    # stephen-boltzman constant
  σ * (T + 273.15)^4
end

function cal_Rln(emiss::FT, T::FT) where {FT<:Real}
  emiss * blackbody(T)
end

# kPa deg-1
function cal_slope(Ta::FT)::FT where {FT<:Real}
  2503.0 / pow((Ta + 237.3), 2) * exp(17.27 * Ta / (Ta + 237.3))
end

# kPa
function cal_es(Ta::FT)::FT where {FT<:Real} 
  0.61078 * exp(17.3 * Ta / (237.3 + Ta))
end

function cal_ea(Ta::FT, RH::FT)::FT where {FT<:Real}
  cal_es(Ta) * RH / 100
end

function cal_lambda(Ta::FT)::FT where {FT<:Real} 
  (2.501 - 0.00237 * Ta) * 1000000
end

function ea2q(ea::FT, Pa::FT=101.35)::FT where {FT<:Real} 
  0.622 * ea / (Pa - 0.378 * ea)
end

function RH2q(Ta::FT, RH::FT)::FT where {FT<:Real}
  es = cal_es(Ta)
  ea = es * RH / 100
  ea2q(ea)
end

"""
# Arguments
- `q`  : specific humidity, g / kg
- `tem`: air temperature, ℃
"""
function q2RH(q::FT, tem::FT)::FT where {FT<:Real}
  # Vapour pressure in mbar
  ea = 0.46 * q * (tem + 273.16) / 100
  es = 6.1078 * exp((17.269 * tem) / (237.3 + tem))
  clamp(ea / es * 100, 0.0, 100.0)
end


function cal_cp(q::FT)::FT where {FT<:Real} 
  1004.65 * (1 + 0.84 * q)
end

function cal_cp(Ta::FT, RH::FT)::FT where {FT<:Real}
  q = RH2q(Ta, RH)
  cal_cp(q)
end

export pow, cal_ea, cal_es,
  cal_slope, cal_lambda, cal_cp, ea2q, RH2q,
  cal_Rln,
  blackbody
