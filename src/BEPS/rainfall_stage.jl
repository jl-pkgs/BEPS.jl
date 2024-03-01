function rainfall_stage1_jl(Tair::Float64, prcp::Float64, 
  pre_mass_water_o::Float64, pre_mass_water_u::Float64,
  lai_o::Float64, lai_u::Float64, clumping::Float64)

  # Ta > 0, otherwise it is snow fall
  Tair <= 0.0 && (prcp = 0.0)
  prcp_o = prcp
  prcp_g = 0.0

  # overstorey
  mass_water_o = pre_mass_water_o + prcp_o * kstep * ρ_w * (1 - exp(-lai_o * clumping))
  massMax_water_o = 0.1 * lai_o

  mass_water_o = clamp(mass_water_o, 0, massMax_water_o)

  massStep_water_o = mass_water_o - pre_mass_water_o
  massStep_water_o = max(0.0, massStep_water_o)

  perc_water_o = mass_water_o / massMax_water_o
  perc_water_o = min(1.0, perc_water_o)

  # understorey
  prcp_u = prcp_o - massStep_water_o / ρ_w / kstep
  mass_water_u = pre_mass_water_u + prcp_u * kstep * ρ_w * (1 - exp(-lai_u * clumping))
  massMax_water_u = 0.1 * lai_u

  mass_water_u = clamp(mass_water_u, 0, massMax_water_u)

  massStep_water_u = mass_water_u - pre_mass_water_u
  massStep_water_u = max(0.0, massStep_water_u)

  perc_water_u = mass_water_u / massMax_water_u
  perc_water_u = min(1, perc_water_u)

  # ground
  prcp_g = prcp_u - massStep_water_u / ρ_w / kstep

  return mass_water_o, mass_water_u, perc_water_o, perc_water_u, prcp_g
end


function rainfall_stage2_jl(evapo_water_o::Float64, evapo_water_u::Float64,
  mass_water_o::Ref{Float64}, mass_water_u::Ref{Float64})

  # kstep  # 6min or 360s per step
  mass_water_o[] = max(mass_water_o[] - evapo_water_o * kstep, 0.0)
  mass_water_u[] = max(mass_water_u[] - evapo_water_u * kstep, 0.0)
end
