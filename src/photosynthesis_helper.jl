function photosynthesis(Tc_old::Leaf, R::Leaf, Ci_old::Leaf, leleaf::Leaf,
  Ta::Cdouble, ea::Cdouble, f_soilwater::Cdouble,
  g0_h2o::Cdouble, g1_h2o::Cdouble,
  Gb_o::Cdouble, Gb_u::Cdouble, Vcmax_sunlit::Cdouble, Vcmax_shaded::Cdouble,
  # output
  Gs_new::Leaf, Ac::Leaf, Ci_new::Leaf; version="c")

  if version == "c"
    fun = photosynthesis_c
  elseif version == "julia"
    fun = photosynthesis_jl
  end

  Gs_new.o_sunlit, Ac.o_sunlit, Ci_new.o_sunlit =
    fun(Tc_old.o_sunlit, R.o_sunlit, ea, Gb_o, Vcmax_sunlit, f_soilwater, g0_h2o, g1_h2o,
      Ci_old.o_sunlit,
      Ta, leleaf.o_sunlit)

  Gs_new.o_shaded, Ac.o_shaded, Ci_new.o_shaded =
    fun(Tc_old.o_shaded, R.o_shaded, ea, Gb_o, Vcmax_shaded, f_soilwater, g0_h2o, g1_h2o,
      Ci_old.o_shaded,
      Ta, leleaf.o_shaded)

  Gs_new.u_sunlit, Ac.u_sunlit, Ci_new.u_sunlit =
    fun(Tc_old.u_sunlit, R.u_sunlit, ea, Gb_u, Vcmax_sunlit, f_soilwater, g0_h2o, g1_h2o,
      Ci_old.u_sunlit,
      Ta, leleaf.u_sunlit)

  Gs_new.u_shaded, Ac.u_shaded, Ci_new.u_shaded =
    fun(Tc_old.u_shaded, R.u_shaded, ea, Gb_u, Vcmax_shaded, f_soilwater, g0_h2o, g1_h2o,
      Ci_old.u_shaded,
      Ta, leleaf.u_shaded)
end


"""
- `ea` : [kPa]
- `Ta` : [°C]
- `ρₐ` : air density [kg m-3]
"""
function cal_rho_a(Ta, ea)
  TK = Ta + 273.13
  ρₐ = ea * 2.165 / TK   # absolute humidity, [kg m-3]
  return ρₐ
end

# atm = 1.013 # 1 atm, [1013.25 hPa] -> 1.013 bar
const pstat273 = 0.022624 / (273.16 * 1.013)

umol_m(gs_mol::T, TK::T) where {T<:Real} = gs_mol * TK * pstat273
m_umol(gs::T, TK::T) where {T<:Real} = gs / (TK * pstat273)

"""
Ball-Berry stomatal conductance model

# Arguments
- `A`  : net photosynthesis rate (mol m-2 s-1)
- `Cs` : CO2 concentration at leaf surface (ppm)
- `RH` : relative humidity at leaf surface (0-1)
- `g0` : the minimum stomatal conductance (umol m-2 s-1)

# Returns
- `gs` : stomatal conductance of co2 (umol m-2 s-1)
"""
function stomatal_conductance(A, Cs, RH, g0, g1, β_soil=1.0)
  gs = g0 + β_soil * g1 * A * RH / Cs
  return gs
end


function sort3(a::T, b::T, c::T) where {T<:Real}
  if a > b
    a, b = b, a
  end
  if b > c
    b, c = c, b
  end
  if a > b
    a, b = b, a
  end
  return a, b, c
end

function findroot(root1::Float64, root2::Float64, root3::Float64)::Float64
  A::Float64 = 0.0
  root_min, root_mid, root_max = sort3(root1, root2, root3)
  # find out where roots plop down relative to the x-y axis
  if root_min > 0 && root_mid > 0 && root_max > 0
    A = root_min
  end
  if root_min < 0 && root_mid < 0 && root_max > 0
    A = root_max
  end
  if root_min < 0 && root_mid > 0 && root_max > 0
    A = root_mid
  end
  A
end


@fastmath function SFC_VPD(T_leaf_K::T, LE::T, λ::T, rᵥ::T, ρₐ::T)::T where {T<:Real}
  es = ES(T_leaf_K)          # mb
  ρᵥ = (LE / λ) * rᵥ + ρₐ    # kg m-3
  e = ρᵥ * T_leaf_K / 0.2165 # mb
  vpd = es - e               # mb
  rh = 1.0 - vpd / es        # 0 to 1.0
  return rh
end

@fastmath function TEMP_FUNC(rate::T, eact::T, tprime::T, tref::T, t_lk::T)::T where {T<:Real}
  rate * exp(tprime * eact / (tref * rugc * t_lk))
end

function LAMBDA(TK::T)::T where {T<:Real}
  y = 3149000.0 - 2370.0 * TK # J kg-1
  # add heat of fusion for melting ice
  if TK < 273.0
    y += 333.0
  end
  return y
end

# Function to calculate saturation vapor pressure function in mb
@fastmath function ES(t::T)::T where {T<:Real}
  if t > 0.0
    y1::T = 54.8781919 - 6790.4985 / t - 5.02808 * log(t)
    y = exp(y1)
  else
    println("bad es calc")
    y = 0.0
  end
  return y
end

# Maxwell-Boltzmann temperature distribution for photosynthesis
@fastmath function TBOLTZ(rate::T, eakin::T, topt::T, tl::T)::T where {T<:Real}
  hkin::T = 200000.0  # enthalpy term, J mol-1

  dtlopt::T = tl - topt
  prodt::T = rugc * topt * tl
  numm::T = hkin * exp(eakin * dtlopt / prodt)
  denom::T = hkin - eakin * (1.0 - exp(hkin * dtlopt / prodt))

  return rate * numm / denom
end
