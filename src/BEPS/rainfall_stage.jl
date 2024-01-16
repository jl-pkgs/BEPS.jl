function rainfall_stage1_jl(Tair::Float64, prcp::Float64, mass_water_o_last::Float64, mass_water_u_last::Float64,
  lai_o::Float64, lai_u::Float64, clumping::Float64)

  length_step = kstep

  # Ta > 0, otherwise it is snow fall
  Tair <= 0.0 && (prcp = 0.0)

  density_water = 1025.0
  prcp_g = 0.0

  # overstorey
  prcp_o = prcp
  mass_water_o = mass_water_o_last + prcp_o * length_step * density_water * (1 - exp(-lai_o * clumping))
  massMax_water_o = 0.1 * lai_o

  mass_water_o = clamp(mass_water_o, 0, massMax_water_o)

  massStep_water_o = mass_water_o - mass_water_o_last
  massStep_water_o = max(0.0, massStep_water_o)

  percent_water_o = mass_water_o / massMax_water_o
  percent_water_o = min(1.0, percent_water_o)

  # understorey
  prcp_u = prcp_o - massStep_water_o / density_water / length_step
  mass_water_u = mass_water_u_last + prcp_u * length_step * density_water * (1 - exp(-lai_u * clumping))
  massMax_water_u = 0.1 * lai_u

  mass_water_u = clamp(mass_water_u, 0, massMax_water_u)

  massStep_water_u = mass_water_u - mass_water_u_last
  massStep_water_u = max(0.0, massStep_water_u)

  percent_water_u = mass_water_u / massMax_water_u
  percent_water_u = min(1, percent_water_u)

  # ground
  prcp_g = prcp_u - massStep_water_u / density_water / length_step

  return mass_water_o, mass_water_u, percent_water_o, percent_water_u, prcp_g
end


function rainfall_stage2_jl(evapo_water_o::Float64, evapo_water_u::Float64,
  mass_water_o::Ref{Float64}, mass_water_u::Ref{Float64})

  length_step = kstep  # 6min or 360s per step
  mass_water_o[] = max(mass_water_o[] - evapo_water_o * length_step, 0.0)
  mass_water_u[] = max(mass_water_u[] - evapo_water_u * length_step, 0.0)
end
