# - `lc`: [1, 2, 6, 9, 13, 40, -1]
function ReadParamVeg(lc::Int=1)
  param = zeros(Float64, 48)
  param[5] = lc

  gen_data = readGeneralParam()
  v = readVegRaw(lc)

  # Vegetation parameters
  param[3] = v["clumping_index"]
  param[9] = v["LAI_max_o"]
  param[10] = v["LAI_max_u"]
  param[17] = v["z00"]
  param[19] = v["mass_overstory"]
  param[21] = v["mass_understory"]
  param[22] = v["root_depth"]
  param[23] = v["albedo_canopy_vis"]
  param[24] = v["albedo_c_nir_observed"]
  param[28] = v["r_root_decay"]
  param[29] = v["minimum_stomatal_resistance"]
  param[30] = v["z_canopy_o"]
  param[31] = v["z_canopy_u"]
  param[34] = v["g1_w"]
  param[37] = v["VCmax25"]
  param[40] = v["leaf_resp_co"] / RTIMES
  param[41] = v["stem_resp_co"] / RTIMES
  param[42] = v["root_resp_co"] / RTIMES
  param[45] = v["fine_root_resp_co"] / RTIMES
  param[47] = v["N_leaf"]
  param[48] = v["slope_Vc"]

  # General parameters
  param[7] = gen_data["light_compensate_point"]
  param[8] = gen_data["light_saturation_point"]
  param[11] = gen_data["LAI_min_overstory"]
  param[12] = gen_data["LAI_min_understory"]
  param[13] = gen_data["albedo_new_snow"]
  param[14] = gen_data["albedo_new_snow_vis"]
  param[15] = gen_data["albedo_new_snow_nir"]
  param[16] = gen_data["density_now_snow"]
  param[18] = gen_data["specific_heat_overstory"]
  param[20] = gen_data["specific_heat_understory"]
  param[25] = gen_data["albedo_saturated_soil"]
  param[26] = gen_data["albedo_dry_soil"]
  param[27] = gen_data["r_drainage"]
  param[32] = gen_data["the_height_to_measure_wind_speed"]
  param[33] = gen_data["the_depth_litter"]
  param[35] = gen_data["intercept_for_H2O_ball_berry"]
  param[36] = gen_data["intercept_for_C2O_ball_berry"]
  param[38] = 2.39 * param[37] - 14.2 # Jmax at 25 C
  param[39] = gen_data["mb_coefficient"]
  param[46] = gen_data["Q10"]

  return param
end


export ReadParamVeg
