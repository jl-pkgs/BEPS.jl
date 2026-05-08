# 附属于`aerodynamic_conductance.jl`而生
@with_kw mutable struct AeroConsts{T}
  ustar::T = 0.0         # friction velocity at reference height [m s-1]
  coef_L::T = 0.0        # coefficient in Monin-Obukhov length: L = coef_L * SH_o_p
  rb_o::T = 200.0        # overstory leaf boundary layer resistance [s m-1]
  rb_u::T = 200.0        # understory leaf boundary layer resistance [s m-1]
  gamma_u::T = 1.0       # attenuation factor for understory wind profile [-]
  exp_u::T = 1.0         # exp(gamma_u * (1 - h_u / h_o)) for ra_u
  exp_g_u::T = 1.0       # exp(4 * (1 - h_u / h_o)) for ground-air resistance
  lai_o_pow075::T = 1.0  # cached lai_o^0.75 term used in gamma_u
  Le_cuberoot::T = 1.0   # cached (lai_o * clumping)^(1/3) for overstory attenuation
end

pow_075(x::T) where {T<:Real} = x^0.75
cuberoot(x::T) where {T<:Real} = x^(1.0 / 3.0)

exp_u_terms(canopy_height_o::T, canopy_height_u::T, gamma_u::T) where {T<:Real} =
  exp(gamma_u * (1 - canopy_height_u / canopy_height_o))

exp_g_terms(canopy_height_o::T, canopy_height_u::T) where {T<:Real} =
  exp(4.0 * (1 - canopy_height_u / canopy_height_o))

"""
    aero_exp_terms(canopy_height_o, canopy_height_u, z_wind, clumping, Tair, wind_sp, lai_o)

Precompute the aerodynamic terms that are independent of the canopy-energy
iteration state (e.g., `SH_o_p` and `L`-dependent stability correction).

Workflow:
1. Compute geometric/wind primitives (`d`, `z0`, `ustar`) from canopy height and wind forcing.
2. Build `coef_L`, which linearly maps sensible heat to Monin-Obukhov length:
   `L = coef_L * SH_o_p`.
3. Precompute expensive nonlinear terms (`lai_o^0.75`, `(lai_o*clumping)^(1/3)`,
   and several `exp(...)` factors) used repeatedly in resistance formulas.
4. Estimate overstory/understory boundary-layer resistances (`rb_o`, `rb_u`)
   from wind profile and Reynolds/Nusselt relationships.

The returned tuple is designed to be cached in `AeroConsts` and reused across
sub-iterations, minimizing repeated `^`/`exp` work while keeping the runtime
path in `aerodynamic_conductance_jl` compact.
"""
function aero_exp_terms(canopy_height_o::T, canopy_height_u::T, z_wind::T, clumping::T,
  Tair::T, wind_sp::T, lai_o::T) where {T<:Real}
  k = 0.4
  cp = 1010.0
  density_air = 1.225
  g = 9.8

  lai_o_pow075 = pow_075(lai_o)
  gamma_u = 0.1 + lai_o_pow075
  exp_u = exp_u_terms(canopy_height_o, canopy_height_u, gamma_u)
  exp_g_u = exp_g_terms(canopy_height_o, canopy_height_u)
  Le = lai_o * clumping
  Le_cuberoot = cuberoot(Le)

  if !(isfinite(wind_sp) && wind_sp > 0)
    return zero(T), zero(T), T(200.0), T(200.0), lai_o_pow075, Le_cuberoot, gamma_u, exp_u, exp_g_u
  end

  # displacement height (m)
  d = 0.8 * canopy_height_o
  # roughness length (m)
  z0 = 0.08 * canopy_height_o
  log_zh_z0 = log((z_wind - d) / z0)
  ustar = wind_sp * k / log_zh_z0 # friction velocity (m/s)
  coef_L = -(k * g) / (density_air * cp * (Tair + 273.3) * ustar^3)

  nu_lower = (13.3 + Tair * 0.07) / 1000000
  alfaw = (18.9 + Tair * 0.07) / 1000000
  uh = 1.1 * ustar / k

  gamma_o = (0.167 + 0.179 * uh) * Le_cuberoot
  ud = uh * exp(-gamma_o * (1 - d / canopy_height_o))
  Nu_o = cal_Nu(ud, nu_lower)
  rb_o = min(40.0, 0.5 * 0.1 / (alfaw * Nu_o))

  un_d = uh * exp(-gamma_u * (1 - canopy_height_u * 0.8 / canopy_height_o))
  Nu_u = cal_Nu(un_d, nu_lower)
  rb_u = min(40.0, 0.5 * 0.1 / (alfaw * Nu_u))

  return ustar, coef_L, rb_o, rb_u, lai_o_pow075, Le_cuberoot, gamma_u, exp_u, exp_g_u
end


function AeroConsts!(ac::AeroConsts{T},
  canopy_height_o::T, canopy_height_u::T, z_wind::T, clumping::T,
  Tair::T, wind_sp::T, lai_o::T) where {T<:Real}

  ustar, coef_L, rb_o, rb_u, lai_o_pow075, Le_cuberoot, gamma_u, exp_u, exp_g_u =
    aero_exp_terms(canopy_height_o, canopy_height_u, z_wind, clumping, Tair, wind_sp, lai_o)

  @pack! ac = ustar, coef_L, rb_o, rb_u, lai_o_pow075, Le_cuberoot, gamma_u, exp_u, exp_g_u
  nothing
end



"""
    ra_updateH(canopy_height_o, canopy_height_u, zh, log_zh_z0, inv_k_ustar, ustar,
      SH_o_p, n, coef_L, gamma_u, exp_u, exp_g_u)

Update aerodynamic resistances (`ra_o`, `ra_u`, `ra_g`) for current sensible
heat flux `SH_o_p`. This function isolates the iteration-dependent stability
part so the main aerodynamic function stays easier to read.
"""
function ra_updateH(SH_o_p::FT, z_wind, canopy_height_o::FT, canopy_height_u::FT,
  ustar::FT, coef_L::FT, gamma_u::FT, exp_u::FT, exp_g_u::FT)

  if !(isfinite(ustar) && ustar > 0 && isfinite(coef_L))
    return FT(200.0), FT(200.0), FT(600.0)
  end

  n = FT(5.0)
  d = 0.8 * canopy_height_o
  z0 = 0.08 * canopy_height_o
  zh = z_wind - d
  log_zh_z0 = log(zh / z0)

  k = 0.4
  inv_k_ustar = 1.0 / (k * ustar)

  L::FT = coef_L * SH_o_p
  L = max(-2.0, L)

  ra_o::FT = inv_k_ustar * (log_zh_z0 + (n * zh * L))
  ra_o = clamp(ra_o, 2, 100)

  if L > 0
    ψ = 1 + 5 * zh * L
  else
    ψ = (1 - 16 * zh * L)^(-0.5)
  end
  ψ = min(10.0, ψ)

  kh_o = 0.41 * ustar * (canopy_height_o - canopy_height_o * 0.8) / ψ
  ra_u = canopy_height_o / (gamma_u * kh_o) * (exp_u - 1)

  # kh_u = kh_o * exp(-4 * (1 - canopy_height_u / canopy_height_o))
  ra_g = canopy_height_o / (4.0 * kh_o) * (exp(4.0) - exp_g_u)
  ra_g = ra_g + ra_u + ra_o
  ra_g = max(120, ra_g)
  return ra_o, ra_u, ra_g
end
