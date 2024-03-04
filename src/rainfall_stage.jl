function rainfall_stage1_jl(Tair::Float64, prcp::Float64,
  perc_water::Layer2{Float64},
  mass_water::Layer2{Float64}, 
  mass_water_pre::Layer2{Float64}, 
  lai_o::Float64, lai_u::Float64, clumping::Float64)

  # Ta > 0, otherwise it is snow fall
  Tair <= 0.0 && (prcp = 0.0)
  prcp_o = prcp
  prcp_g = 0.0

  # overstorey
  mass_water.o = mass_water_pre.o + prcp_o * kstep * ρ_w * (1 - exp(-lai_o * clumping))
  massMax_water_o = 0.1 * lai_o

  mass_water.o = clamp(mass_water.o, 0, massMax_water_o)
  massStep_water_o = max(mass_water.o - mass_water_pre.o, 0.0)

  perc_water.o = min(mass_water.o / massMax_water_o, 1.0)

  # understorey
  prcp_u = prcp_o - massStep_water_o / ρ_w / kstep
  mass_water.u = mass_water_pre.u + prcp_u * kstep * ρ_w * (1 - exp(-lai_u * clumping))
  massMax_water_u = 0.1 * lai_u

  mass_water.u = clamp(mass_water.u, 0, massMax_water_u)
  massStep_water_u = max(mass_water.u - mass_water_pre.u, 0.0)

  perc_water.u = min(mass_water.u / massMax_water_u, 1.0)
  prcp_g = prcp_u - massStep_water_u / ρ_w / kstep

  return prcp_g
end


function rainfall_stage2_jl(evapo_water_o::Float64, evapo_water_u::Float64,
  mass_water::Layer2{Float64})

  # kstep  # 6min or 360s per step
  mass_water.o = max(mass_water.o - evapo_water_o * kstep, 0.0)
  mass_water.u = max(mass_water.u - evapo_water_u * kstep, 0.0)
end
