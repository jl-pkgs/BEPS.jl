## soil evap
perc_snow_g = var.Xcs_g[kkk]
depth_water = depth_water[]
depth_snow = depth_snow[]
mass_water_g = ρ_w * depth_water[]
mass_snow_g = var.Wcs_g[kkk]
m = SurfaceMass(; perc_snow_g, depth_water, depth_snow, mass_water_g, mass_snow_g)

var.Evap_soil[kkk], var.Evap_SW[kkk], var.Evap_SS[kkk] =
  evaporation_soil_jl(temp_grd, var.Ts0[kkk-1], rh_air, radiation_g, Gheat_g,
    m,
    # Ref(var.Xcs_g, kkk), depth_water, depth_snow, mass_water_g, Ref(var.Wcs_g, kkk), # Ref
    var.rho_snow[kkk], soilp.θ_prev[1], soilp.θ_sat[1])
# Ref(var.Evap_soil, kkk), Ref(var.Evap_SW, kkk), Ref(var.Evap_SS, kkk)

var.Xcs_g[kkk] = m.perc_snow_g
depth_water[] = m.depth_water
depth_snow[] = m.depth_snow
mass_water_g = m.mass_water_g
var.Wcs_g[kkk] = m.mass_snow_g


# var.Evap_soil[kkk], var.Evap_SW[kkk], var.Evap_SS[kkk] =
#   evaporation_soil_jl(temp_grd, var.Ts0[kkk-1], rh_air, radiation_g, Gheat_g,
#     Ref(var.Xcs_g, kkk), depth_water, depth_snow, mass_water_g, Ref(var.Wcs_g, kkk), # Ref
#     var.rho_snow[kkk], soilp.θ_prev[1], soilp.θ_sat[1])
