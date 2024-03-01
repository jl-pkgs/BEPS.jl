## soil evap
perc_snow_g = var.Xg_snow[kkk]
depth_water = Zp[]
depth_snow = Zsp[]
mass_water_g = ρ_w * Zp[]
mass_snow_g = var.Wg_snow[kkk]
m = SurfaceMass(; perc_snow_g, depth_water, depth_snow, mass_water_g, mass_snow_g)

var.Evap_soil[kkk], var.Evap_SW[kkk], var.Evap_SS[kkk] =
  evaporation_soil_jl(temp_grd, var.Ts0[kkk-1], rh_air, radiation_g, Gheat_g,
    m,
    # Ref(var.Xg_snow, kkk), Zp, Zsp, mass_water_g, Ref(var.Wg_snow, kkk), # Ref
    var.rho_snow[kkk], soilp.θ_prev[1], soilp.θ_sat[1])
# Ref(var.Evap_soil, kkk), Ref(var.Evap_SW, kkk), Ref(var.Evap_SS, kkk)

var.Xg_snow[kkk] = m.perc_snow_g
Zp[] = m.depth_water
Zsp[] = m.depth_snow
mass_water_g = m.mass_water_g
var.Wg_snow[kkk] = m.mass_snow_g


# var.Evap_soil[kkk], var.Evap_SW[kkk], var.Evap_SS[kkk] =
#   evaporation_soil_jl(temp_grd, var.Ts0[kkk-1], rh_air, radiation_g, Gheat_g,
#     Ref(var.Xg_snow, kkk), Zp, Zsp, mass_water_g, Ref(var.Wg_snow, kkk), # Ref
#     var.rho_snow[kkk], soilp.θ_prev[1], soilp.θ_sat[1])
