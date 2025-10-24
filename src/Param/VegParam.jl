export VegParam, theta2par, theta2par!, par2theta, par2theta!


@with_kw mutable struct VegParam{FT<:AbstractFloat}
  # lc::Int = 1
  # Ω0::FT = 0.7             # // clumping_index
  LAI_max_o::FT = 4.5     # ? 应该不重要
  LAI_max_u::FT = 2.4     

  α_canopy_vis::FT = 0.055 # 0.04
  α_canopy_nir::FT = 0.300 # 0.25
  α_soil_sat::FT = 0.10    # albedo of saturated/dry soil, `rainfall1`
  α_soil_dry::FT = 0.35    # the albedo of dry soil

  # r_drainage::FT = 0.5     # ? 产流比例
  # r_root_decay::FT = 0.97  # ? decay_rate_of_root_distribution

  z_canopy_o::FT = 1.0      # [m]
  z_canopy_u::FT = 0.2      # [m]
  z_wind::FT = 2            # [m]

  g1_w::FT = 8             # Ball-Berry, slope coefficient
  g0_w::FT = 0.0175        # Ball-Berry, intercept_for_H2O

  VCmax25::FT = 89.45       # maximum capacity of Rubisco at 25℃
  # Jmax25::FT = 2.39 * 57.7 - 14.2 # 

  # # coefficient reflecting the sensitivity of stomata to VPD/moderately N-stressed plants
  N_leaf::FT = 1.74 + 0.71   # leaf Nitrogen content, mean value + 1 SD [g/m2]
  slope_Vc::FT = 33.79 / 57.7
end


function theta2par(theta::Vector{FT}) where {FT<:Real}
  par = VegParam{FT}()
  theta2par!(par, theta)
end

function theta2par!(par::VegParam{FT}, theta::Vector{FT}) where {FT<:Real}
  # par.Ω0 = theta[3]
  # par.lc = theta[5]
  par.LAI_max_o = theta[9]
  par.LAI_max_u = theta[10]
  par.α_canopy_vis = theta[23]
  par.α_canopy_nir = theta[24]
  par.α_soil_sat = theta[25]
  par.α_soil_dry = theta[26]
  # par.r_drainage = theta[27]
  # par.r_root_decay = theta[28]
  par.z_canopy_o = theta[30]
  par.z_canopy_u = theta[31]
  par.z_wind = theta[32]
  par.g1_w = theta[34]
  par.g0_w = theta[35]
  par.VCmax25 = theta[37]
  par.N_leaf = theta[47]
  par.slope_Vc = theta[48]
  return par
end

function par2theta(par::VegParam{FT}; kw...) where {FT<:Real}
  theta = zeros(FT, 48)
  par2theta!(theta, par; kw...)
end

function par2theta!(theta::Vector{FT}, par::VegParam{FT}; clumping=0.7, VegType::Int=1)
  theta[3] = clumping
  theta[5] = VegType
  theta[9] = par.LAI_max_o
  theta[10] = par.LAI_max_u
  theta[23] = par.α_canopy_vis
  theta[24] = par.α_canopy_nir
  theta[25] = par.α_soil_sat
  theta[26] = par.α_soil_dry
  # theta[27] = par.r_drainage
  # theta[28] = par.r_root_decay
  theta[30] = par.z_canopy_o
  theta[31] = par.z_canopy_u
  theta[32] = par.z_wind
  theta[34] = par.g1_w
  theta[35] = par.g0_w
  theta[37] = par.VCmax25
  theta[47] = par.N_leaf
  theta[48] = par.slope_Vc
  return theta
end

# inds = [
#   3, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
#   26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 45, 46, 47, 48
# ]

## Removed Parameters
# I_com::FT = 100          # light_compensate_point 
# I_sat::FT = 1000         # light_saturation_point, []

# mb::FT = 28.0
# z_litter::FT = 0.05    # [m]
# z_root::FT = 0.8       # [m]
# rs_min::FT = 200       # minimum_stomatal_resistance, [s m-1]

# LAI_min_o::FT = 0.01
# LAI_min_u::FT = 0.01

# z00::FT = 1.53           # roughness length for heat
# cp_o::FT = 2700          # specific_heat, [J kg-1 K-1]
# cp_u::FT = 2700          # specific_heat, [J kg-1 K-1]
# mass_o::FT = 40          # [kg m-2]
# mass_u::FT = 10          # [kg m-2]
# ρ_new_snow::FT = 100.0   # [kg m-3]

# albedo
# α_snow::FT = 0.87        # new snow
# α_snow_vis::FT = 0.94
# α_snow_nir::FT = 0.8

# leaf_resp_co::FT = 0.015 / RTIMES      # kg C-1 d-1 kg-1
# stem_resp_co::FT = 0.0035 / RTIMES     # kg C-1 d-1 kg-1
# root_resp_co::FT = 0.0025 / RTIMES     # kg C-1 d-1 kg-1
# fine_root_resp_co::FT = 0.003 / RTIMES # kg C-1 d-1 kg-1
# Q10::FT = 2.3                          # constant for exp. resp.
