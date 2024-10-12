include("photosynthesis_helper.jl")

"""
A = Ag - Rd, net photosynthesis is the difference between gross photosynthesis
and dark respiration. Note photorespiration is already factored into Ag.

Gs from Ball-Berry is for water vapor.  It must be divided by the ratio of the
molecular diffusivities to be valid for A.

Forests are hypostomatous. Hence, we don't divide the total resistance by 2
since transfer is going on only one side of a leaf.

# Arguments
- `ea`     : [kPa]
- `gb_w`   : leaf laminar boundary layer conductance to H2O, [s m-1]
- `Vcmax25`: the maximum rate of carboxylation of Rubisco at 25 deg C, [umol m-2
  s-1]

- `cii`   : intercellular CO2 concentration (ppm)
- `g0_h2o`: the minimum stomatal conductance to H2O,      [umol m-2 s-1]
- `g1_h2o`: the slope of the stomatal conductance to H2O, [unitless]

- `LH_leaf`: latent heat of vaporization of water at the leaf temperature, [W
  m-2]

# Intermediate variables
- `gb_c`     : leaf laminar boundary layer condunctance to CO2 (mol m-2 s-1)
- `RH_leaf`    : relative humidity at leaf surface (0-1)
- `gs_co2_mole`: stomatal conductance to CO2 (mol m-2 s-1)
- `gs_h2o_mole`: stomatal conductance to h2o (mol m-2 s-1)
- `cs`         : CO2 concentration at leaf surface (ppm)
- `b_co2`      : the intercept term in BWB model (mol CO2 m-2 s-1): b_h2o/1.6
- `m_co2`      : the slope in BWB model: m_h2o/1.6
- `Γ`          : CO2 compensation point (ppm)
- `jmopt`      : the maximum potential electron transport rate at 25 deg C (umol
  m-2 s-1)
- `Jmax`       : the maximum potential electron transport rate (umol m-2 s-1)
- `Vcmax`      : the maximum velocities of carboxylation of Rubisco (umol m-2
  s-1)
- `km_co2`     : Michaelis-Menten constant for CO2 (µmol mol-1)
- `km_o2`      : Michaelis-Menten constant for O2 (mmol mol-1)
- `tau`        : the specifity of Rubisco for CO2 compared with O2
- `Jₓ`         : the flux of electrons through the thylakoid membrane (umol m-2
  s-1)
"""
@fastmath function photosynthesis_jl(T_leaf_p::Cdouble, Rsn_leaf::Cdouble, ea::Cdouble,
  gb_w::Cdouble, Vcmax25::Cdouble,
  β_soil::Cdouble, g0_w::Cdouble, g1_w::Cdouble,
  cii::Cdouble,
  T_leaf::Cdouble, LH_leaf::Cdouble)

  ca::Cdouble = CO2_air                 # atmospheric co2 concentration (ppm)
  PPFD::Cdouble = 4.55 * 0.5 * Rsn_leaf # incident photosynthetic photon flux density (PPFD) umol m-2 s-1
  (2PPFD < 1) && (PPFD = 0.0)
  T_leaf_K = T_leaf + 273.13

  λ = LAMBDA(T_leaf_p) # [J kg-1], ~2.5MJ kg-1
  rᵥ = 1.0 / gb_w

  ρₐ = cal_rho_a(T_leaf, ea) # [kg m-3]
  gb_c_mol = m_umol(gb_w / 1.6, T_leaf_K) # [s m-1] -> [umol m-2 s-1]

  g0_c = g0_w / 1.6
  g1_c = g1_w / 1.6

  RH_leaf = SFC_VPD(T_leaf_K, LH_leaf, λ, rᵥ, ρₐ)
  tprime25 = T_leaf_K - TK25  # temperature difference

  Kc = TEMP_FUNC(kc25, ekc, tprime25, TK25, T_leaf_K)
  Ko = TEMP_FUNC(ko25, eko, tprime25, TK25, T_leaf_K)
  tau = TEMP_FUNC(tau25, ektau, tprime25, TK25, T_leaf_K)

  K = Kc * (1.0 + o2 / Ko)
  Γ = 0.5 * o2 / tau * 1000.0  # umol mol-1
  Rd25 = Vcmax25 * 0.004657    # leaf dark respiration (umol m-2 s-1)

  # Bin Chen: Reduce respiration by 40% in light according to Amthor
  (2PPFD > 10) && (Rd25 *= 0.4)
  Rd = TEMP_FUNC(Rd25, erd, tprime25, TK25, T_leaf_K)

  #	jmopt = 29.1 + 1.64*Vcmax25; Chen 1999, Eq. 7
  jmopt = 2.39 * Vcmax25 - 14.2

  Jmax = TBOLTZ(jmopt, ejm, toptjm, T_leaf_K)    # Apply temperature correction to JMAX
  Vcmax = TBOLTZ(Vcmax25, evc, toptvc, T_leaf_K)  # Apply temperature correction to vcmax

  #  * Test for the minimum of Wc and Wj.  Both have the form: `W = (a ci - a d)/(e ci + b)`
  #  * after the minimum is chosen set a, b, e and d for the cubic solution.
  #  * estimate of J according to Farquhar and von Cammerer (1981)
  # /*if (jmax > 0) Jₓ = qalpha * iphoton / sqrt(1. +(qalpha2 * iphoton * iphoton / (jmax * jmax)));
  Jₓ = Jmax * PPFD / (PPFD + 2.1 * Jmax) # chen1999, eq.6, J photon from Harley

  # initial guess of intercellular CO2 concentration to estimate Wc and Wj:
  wj = Jₓ * (cii - Γ) / (4.0 * cii + 8.0 * Γ)
  wc = Vcmax * (cii - Γ) / (cii + K)

  if (wj < wc)
    # A_guess = wj      # Harley and Farquhar type model for Wj
    a = Jₓ
    b = 8.0 * Γ
    e = 4.0
  else
    # A_guess = wc
    a = Vcmax # x1
    b = K
    e = 1.0
  end
  d = Γ
  An = 0.0

  if !(wj <= Rd || wc <= Rd)
    # g_s =  g0 + g1 * RH_leaf * β_soil * A_g
    # g_s =  g0 + _c * A_g, (c = g1 * RH_leaf * β_soil)
    c = g1_c * RH_leaf * β_soil
    α = 1.0 + (g0_c / gb_c_mol) - c
    β = ca * (gb_c_mol * c - 2.0 * g0_c - gb_c_mol)
    γ = ca * ca * gb_c_mol * g0_c
    θ = gb_c_mol * c - g0_c

    An = solve_cubic(α, β, γ, θ, Rd, a, b, d, e, ca)
    # Sucrose limitation of photosynthesis, as suggested by Collatz.  `Js=Vmax/2`
    # net photosynthesis rate limited by sucrose synthesis (umol m-2 s-1)
    j_sucrose = Vcmax / 2.0 - Rd
    An = min(An, j_sucrose)
  end

  if An <= 0.0
    An = solve_quad(ca, gb_c_mol, g0_c, a, b, d, e, Rd)
  end
  An = max(0.0, An)
  cs = ca - An / gb_c_mol

  gs_w_mol = (β_soil * g1_w * RH_leaf * An / cs) + g0_w  # mol m-2 s-1
  gs_c_mol = gs_w_mol / 1.6

  ci = cs - An / gs_c_mol
  gs_w = umol_m(gs_w_mol, T_leaf_K)  # s m-1
  return gs_w, An, ci
end


"""
If `wj` or `wc` are less than Rd then A would probably be less than 0. This
would yield a negative stomatal conductance. 

In this case, assume `gs` equals the cuticular value `g0`. This assumptions
yields a quadratic rather than cubic solution for A.

If `A < 0`, set stomatal conductance to cuticle value. A quadratic solution of A
is derived if gs=b, but a cubic form occur if gs = ax + b.  Use quadratic case
when `A<=0`.

```math
c_s = c_a - A * 1/ g_b
c_i = c_a - A * (1/g_b + 1/g_s)
g_s = g_0
Ag = a(c_i - Gamma) / (e c_i + b})
A = Ag - Rd
```
"""
function solve_quad(ca::T, gb_c_mol::T, g0_c::T, a::T, b::T, d::T, e::T, Rd::T) where {T<:Real}
  D = 1 / gb_c_mol + 1 / g0_c
  _a = -D * e
  _b = e * ca + b - e * Rd * D + a * D
  _c = b * Rd - a * ca + a * d + e * Rd * ca

  Δ = _b^2 - 4 * _a * _c
  if Δ >= 0.0
    An = (-_b + sqrt(Δ)) / 2_a
  else
    An = 0.0
  end
  return An
end


"""
Cubic solution: `A^3 + p A^2 + q A + r = 0`. Let `A = x - p / 3`, => `x^3 + ax + b = 0`
Rank roots #1, #2 and #3 according to the minimum, intermediate and maximum value

```math
c_s = c_a - A * 1/ g_b
c_i = c_a - A * (1/g_b + 1/g_s)
g_s = g_0 + g_1 RH A / c_s
Ag = a(c_i - Gamma) / (e c_i + b})
A = Ag - Rd
```
"""
function solve_cubic(α::T, β::T, γ::T, θ::T, Rd::T, a::T, b::T, d::T, e::T, ca::T) where {T<:Real}
  m = e * α
  p = (e * β + b * θ - a * α + e * Rd * α) / m
  q = (e * γ + (b * γ / ca) - a * β + a * d * θ + e * Rd * β + Rd * b * θ) / m
  r = (-a * γ + a * d * γ / ca + e * Rd * γ + Rd * b * γ / ca) / m

  # Use solution from Numerical Recipes from Press
  Q = (p * p - 3.0 * q) / 9.0
  U = (2.0 * p * p * p - 9.0 * p * q + 27.0 * r) / 54.0
  (Q < 0) && return T(0.0)

  r3q = U / sqrt(Q * Q * Q)
  r3q = clamp(r3q, -1, 1) #  by G. Mo
  ψ = acos(r3q)

  root1 = -2sqrt(Q) * cos(ψ / 3.0) - p / 3.0  # real roots
  root2 = -2sqrt(Q) * cos((ψ + 2pi) / 3.0) - p / 3.0
  root3 = -2sqrt(Q) * cos((ψ - 2pi) / 3.0) - p / 3.0
  An = findroot(root1, root2, root3)
  return An
end
