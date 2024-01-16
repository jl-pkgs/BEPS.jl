function Init_Soil_Parameters(landcover::Integer, stxt::Integer, r_root_decay::Float64, p::Soil)
  p.n_layer = 5

  if landcover == 6 || landcover == 9
    p.psi_min = 10.0 # ψ_min
    p.alpha = 1.5
  else
    p.psi_min = 33.0 # ψ_min
    p.alpha = 0.4
  end

  p.d_soil = [0.05, 0.10, 0.20, 0.40, 1.25]

  p.r_root_decay = r_root_decay
  SoilRootFraction(p)

  p.density_soil = [1300.0, 1500.0, 1517.0, 1517.0, 1517.0] # from flux tower.
  p.f_org = [5, 2, 1, 1, 0.3]

  if stxt == 1  # sand
    p.b = [1.7, 1.9, 2.1, 2.3, 2.5]
    p.Ksat = [0.000058, 0.000052, 0.000046, 0.000035, 0.000010]  # saturated hydraulic conductivity
    p.fei = [0.437, 0.437, 0.437, 0.437, 0.437]  # porosity
    p.theta_vfc = [0.09, 0.09, 0.09, 0.09, 0.09]  # field capacity
    p.theta_vwp = [0.03, 0.03, 0.03, 0.03, 0.03]  # wilt point
    p.thermal_cond = [8.6, 8.6, 8.6, 8.6, 8.6]  # thermal conductivity
    p.psi_sat = [0.07, 0.08, 0.09, 0.10, 0.12]  # water potential at sat

  elseif stxt == 2  # loamy sand
    p.b = [2.1, 2.3, 2.5, 2.7, 2.9]
    p.Ksat = [0.000017, 0.000015, 0.000014, 0.000010, 0.000003]  # saturated hydraulic conductivity
    p.fei = [0.437, 0.437, 0.437, 0.437, 0.437]  # porosity
    p.theta_vfc = [0.21, 0.21, 0.21, 0.21, 0.21]  # field capacity
    p.theta_vwp = [0.06, 0.06, 0.06, 0.06, 0.06]  # wilt point
    p.thermal_cond = [8.3, 8.3, 8.3, 8.3, 8.3]  # thermal conductivity
    p.psi_sat = [0.09, 0.10, 0.11, 0.12, 0.14]  # water potential at sat

  elseif stxt == 3  # sandy loam
    p.b = [3.1, 3.3, 3.5, 3.7, 3.9]
    p.Ksat = [0.0000072, 0.00000648, 0.00000576, 0.00000432, 0.00000144]  # saturated hydraulic conductivity
    p.fei = [0.453, 0.453, 0.453, 0.453, 0.453]  # porosity
    p.theta_vfc = [0.21, 0.21, 0.21, 0.21, 0.21]  # field capacity
    p.theta_vwp = [0.10, 0.10, 0.10, 0.10, 0.10]  # wilt point
    p.thermal_cond = [8.0, 8.0, 8.0, 8.0, 8.0]  # thermal conductivity
    p.psi_sat = [0.15, 0.16, 0.17, 0.18, 0.20]  # water potential at sat

  elseif stxt == 4  # loam
    p.b = [4.5, 4.7, 4.9, 5.1, 5.3]
    p.Ksat = [0.0000037, 0.0000033, 0.00000296, 0.00000222, 0.00000074]  # saturated hydraulic conductivity
    p.fei = [0.463, 0.463, 0.463, 0.463, 0.463]  # porosity
    p.theta_vfc = [0.27, 0.27, 0.27, 0.27, 0.27]  # field capacity
    p.theta_vwp = [0.12, 0.12, 0.12, 0.12, 0.12]  # wilt point
    p.thermal_cond = [7.0, 7.0, 7.0, 7.0, 7.0]  # thermal conductivity
    p.psi_sat = [0.11, 0.12, 0.13, 0.14, 0.16]  # water potential at sat

  elseif stxt == 5  # silty loam
    p.b = [4.7, 4.9, 5.1, 5.3, 5.5]
    p.Ksat = [0.0000019, 0.0000017, 0.00000152, 0.00000114, 0.00000038]  # saturated hydraulic conductivity
    p.fei = [0.501, 0.501, 0.501, 0.501, 0.501]  # porosity
    p.theta_vfc = [0.33, 0.33, 0.33, 0.33, 0.33]  # field capacity
    p.theta_vwp = [0.13, 0.13, 0.13, 0.13, 0.13]  # wilt point
    p.thermal_cond = [6.3, 6.3, 6.3, 6.3, 6.3]  # thermal conductivity
    p.psi_sat = [0.21, 0.22, 0.23, 0.24, 0.26]  # water potential at sat

  elseif stxt == 6 # sandy clay loam
    p.b = [4.0, 4.2, 4.4, 4.6, 4.8]
    p.Ksat = [0.0000012, 0.00000108, 0.0000096, 0.0000072, 0.0000024]
    p.fei = fill(0.398, 5)
    p.theta_vfc = fill(0.26, 5)
    p.theta_vwp = fill(0.15, 5)
    p.thermal_cond = fill(7.0, 5)
    p.psi_sat = [0.28, 0.29, 0.30, 0.31, 0.33]

  elseif stxt == 7 # clay loam
    p.b = [5.2, 5.4, 5.6, 5.8, 6.0]
    p.Ksat = [0.00000064, 0.00000058, 0.00000051, 0.00000038, 0.00000013]
    p.fei = fill(0.464, 5)
    p.theta_vfc = fill(0.32, 5)
    p.theta_vwp = fill(0.20, 5)
    p.thermal_cond = [5.8, 5.8, 5.7, 5.8, 5.8]
    p.psi_sat = [0.26, 0.27, 0.28, 0.29, 0.31]

  elseif stxt == 8 # silty clay loam
    p.b = [6.6, 6.8, 7.0, 7.2, 7.4]
    p.Ksat = [0.00000042, 0.00000038, 0.00000034, 0.000000252, 0.000000084]
    p.fei = fill(0.471, 5)
    p.theta_vfc = fill(0.37, 5)
    p.theta_vwp = fill(0.32, 5)
    p.thermal_cond = fill(4.2, 5)
    p.psi_sat = [0.33, 0.34, 0.35, 0.36, 0.38]

  elseif stxt == 9 # sandy clay
    p.b = [6.0, 6.2, 6.4, 6.6, 6.8]
    p.Ksat = [0.00000033, 0.0000003, 0.000000264, 0.000000198, 0.000000066]
    p.fei = fill(0.430, 5)
    p.theta_vfc = fill(0.34, 5)
    p.theta_vwp = fill(0.24, 5)
    p.thermal_cond = fill(6.3, 5)
    p.psi_sat = [0.29, 0.30, 0.31, 0.32, 0.34]

  elseif stxt == 10 # silty clay
    p.b = [7.9, 8.1, 8.3, 8.5, 8.7]
    p.Ksat = [0.00000025, 0.000000225, 0.0000002, 0.00000015, 0.00000005]
    p.fei = fill(0.479, 5)
    p.theta_vfc = fill(0.39, 5)
    p.theta_vwp = fill(0.25, 5)
    p.thermal_cond = fill(4.0, 5)
    p.psi_sat = [0.34, 0.35, 0.36, 0.37, 0.39]

  elseif stxt == 11 # clay
    p.b = [7.6, 7.8, 8.0, 8.2, 8.4]
    p.Ksat = [0.00000017, 0.000000153, 0.000000136, 0.000000102, 0.000000034]
    p.fei = fill(0.475, 5)
    p.theta_vfc = fill(0.40, 5)
    p.theta_vwp = fill(0.27, 5)
    p.thermal_cond = fill(4.4, 5)
    p.psi_sat = [0.37, 0.38, 0.39, 0.40, 0.42]

  else # default
    p.b = [7.6, 7.8, 8.0, 8.2, 8.4]
    p.Ksat = [0.00000017, 0.000000153, 0.000000136, 0.000000102, 0.000000034]
    p.fei = fill(0.475, 5)
    p.theta_vfc = fill(0.40, 5)
    p.theta_vwp = fill(0.27, 5)
    p.thermal_cond = fill(4.4, 5)
    p.psi_sat = [0.37, 0.38, 0.39, 0.40, 0.42]
  end
  p
end


"""
update `f_root`
"""
function SoilRootFraction(soil::Soil)
  cum_depth = zeros(soil.n_layer)

  # For the first layer
  cum_depth[1] = soil.d_soil[1]
  soil.f_root[1] = 1 - (soil.r_root_decay)^(cum_depth[1] * 100)

  # For 2 to n_layer - 1
  for i in 2:(soil.n_layer-1)
    cum_depth[i] = cum_depth[i-1] + soil.d_soil[i]
    soil.f_root[i] = (soil.r_root_decay)^(cum_depth[i-1] * 100) - (soil.r_root_decay)^(cum_depth[i] * 100)
  end

  # For the last layer
  soil.f_root[soil.n_layer] = (soil.r_root_decay)^(cum_depth[soil.n_layer-1] * 100)
end


function Init_Soil_Status(p::Soil, Tsoil::Float64, Tair::Float64, Ms::Float64, snowdepth::Float64)
  d_t = clamp(Tsoil - Tair, -5.0, 5.0)

  p.Zp = 0.0
  p.Zsp = snowdepth
  p.r_rain_g = 0.0

  temp_scale_factors = [0.4, 0.5, 1.0, 1.2, 1.4]
  moisture_scale_factors = [0.8, 1.0, 1.05, 1.10, 1.15]

  for i in 1:5
    p.temp_soil_c[i] = Tair + temp_scale_factors[i] * d_t
    p.temp_soil_p[i] = Tair + temp_scale_factors[i] * d_t
    p.thetam[i] = moisture_scale_factors[i] * Ms
    p.thetam_prev[i] = moisture_scale_factors[i] * Ms
  end

  for i in 1:p.n_layer
    if p.temp_soil_c[i] < -1.0
      p.ice_ratio[i] = 1.0
    elseif p.temp_soil_c[i] > 0
      p.ice_ratio[i] = 0
    else
      p.ice_ratio[i] = -p.temp_soil_c[i] / 1.0
    end
  end
end

function Update_temp_soil_c(p::Soil, value::Cdouble)
  p.temp_soil_c[1] = value
end

function Update_G(p::Soil, value::Cdouble)
  p.G[1] = value
end

export Init_Soil_Parameters, Init_Soil_Status, SoilRootFraction, Update_temp_soil_c, Update_G
