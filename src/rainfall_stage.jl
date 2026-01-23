# - m    : [kg m-2]
# - prcp : [m m-2 s-1]
function water_change(m_water_pre, prcp, lai, Ω)
  mMax_water = 0.1 * lai
  τ = 1 - exp(-lai * Ω)
  m_water_o = m_water_pre + prcp * kstep * ρ_w * τ
  m_water_o = clamp(m_water_o, 0, mMax_water)

  Δm_water_o = max(m_water_o - m_water_pre, 0.0)
  frac_water_o = min(m_water_o / mMax_water, 1.0)
  m_water_o, frac_water_o, Δm_water_o
end

# [kg m-2] -> [m s-1]
kg2m(p) = p / ρ_w / kstep
m2kg(m) = m * kstep * ρ_w

# - m_water: change
function rainfall_stage1_jl(Tair::Float64, prcp::Float64,
  frac_water::Layer2{Float64}, m_water::Layer2{Float64}, m_water_pre::Layer2{Float64},
  lai_o::Float64, lai_u::Float64, Ω::Float64)
  # Ta > 0, otherwise it is snow fall
  Tair <= 0.0 && (prcp = 0.0)
  prcp_o = prcp

  # overstorey
  m_water.o, frac_water.o, Δm_water_o = water_change(m_water_pre.o, prcp_o, lai_o, Ω)
  # understorey
  prcp_u = prcp_o - Δm_water_o / ρ_w / kstep
  m_water.u, frac_water.u, Δm_water_u = water_change(m_water_pre.u, prcp_u, lai_u, Ω)

  prcp_g = prcp_u - Δm_water_u / ρ_w / kstep
  return prcp_g
end

function rainfall_stage2_jl(evapo_water_o::Float64, evapo_water_u::Float64,
  mass_water::Layer2{Float64})
  mass_water.o = max(mass_water.o - evapo_water_o * kstep, 0.0)
  mass_water.u = max(mass_water.u - evapo_water_u * kstep, 0.0)
end
