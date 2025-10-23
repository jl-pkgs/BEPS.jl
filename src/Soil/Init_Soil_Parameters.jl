"""
    Init_Soil_Parameters(landcover::Integer, stxt::Integer, r_root_decay::Float64, p::Soil)

Initialize soil parameters

- `Ksat`         : saturated hydraulic conductivity
- `porosity`     : porosity
- `θ_vfc`        : field capacity
- `θ_vwp`        : wilt point
- `ψ_sat`        : water potential at saturate
- `κ_dry`            : thermal conductivity
"""
function Init_Soil_Parameters(VegType::Integer, SoilType::Integer, r_root_decay::Float64, p::Soil)
  p.n_layer = 5

  if VegType == 6 || VegType == 9 # DBF or EBF, low constaint threshold
    p.ψ_min = 10.0 # ψ_min
    p.alpha = 1.5
  else
    p.ψ_min = 33.0 # ψ_min
    p.alpha = 0.4
  end

  p.dz[1:5] .= [0.05, 0.10, 0.20, 0.40, 1.25] # BEPS V2023
  # z = [0, 5, 15, 25, 35, 45, 55.0] ./ 100
  # z_mid = (z[1:end-1] .+ z[2:end]) ./ 2
  # dz = diff(z)
  # n = length(z)
  # p.dz[1:n-1] = dz

  p.r_root_decay = r_root_decay
  SoilRootFraction(p)

  p.density_soil[1:5] .= [1300.0, 1500.0, 1517.0, 1517.0, 1517.0] # from flux tower.
  p.f_org[1:5] .= [5, 2, 1, 1, 0.3]

  n = 5
  if SoilType == 1  # sand
    b = [1.7, 1.9, 2.1, 2.3, 2.5]
    Ksat = [0.000058, 0.000052, 0.000046, 0.000035, 0.000010]  # 
    porosity = fill(0.437, n) # porosity
    θ_vfc = fill(0.09, n)     # field capacity
    θ_vwp = fill(0.03, n)     # wilt point
    κ_dry = fill(8.6, n) # thermal conductivity
    ψ_sat = [0.07, 0.08, 0.09, 0.10, 0.12]   # water potential at sat

  elseif SoilType == 2  # loamy sand
    b = [2.1, 2.3, 2.5, 2.7, 2.9]
    Ksat = [0.000017, 0.000015, 0.000014, 0.000010, 0.000003] 
    porosity = fill(0.437, n)  # porosity
    θ_vfc = fill(0.21, n)  # field capacity
    θ_vwp = fill(0.06, n)  # wilt point
    κ_dry = fill(8.3, n)  # thermal conductivity
    ψ_sat = [0.09, 0.10, 0.11, 0.12, 0.14]  # water potential at sat

  elseif SoilType == 3  # sandy loam
    b = [3.1, 3.3, 3.5, 3.7, 3.9]
    Ksat = [0.0000072, 0.00000648, 0.00000576, 0.00000432, 0.00000144]
    porosity = fill(0.453, n)  # porosity
    θ_vfc = fill(0.21, n)  # field capacity
    θ_vwp = fill(0.10, n)  # wilt point
    κ_dry = fill(8.0, n)  # thermal conductivity
    ψ_sat = [0.15, 0.16, 0.17, 0.18, 0.20]  # water potential at sat

  elseif SoilType == 4  # loam
    b = [4.5, 4.7, 4.9, 5.1, 5.3]
    Ksat = [0.0000037, 0.0000033, 0.00000296, 0.00000222, 0.00000074] 
    porosity = [0.463, 0.463, 0.463, 0.463, 0.463]  # porosity
    θ_vfc = fill(0.27, n)  # field capacity
    θ_vwp = fill(0.12, n)  # wilt point
    κ_dry = fill(7.0, n)  # thermal conductivity
    ψ_sat = [0.11, 0.12, 0.13, 0.14, 0.16]  # water potential at sat

  elseif SoilType == 5  # silty loam
    b = [4.7, 4.9, 5.1, 5.3, 5.5]
    Ksat = [0.0000019, 0.0000017, 0.00000152, 0.00000114, 0.00000038] 
    porosity = fill(0.501, n)  # porosity
    θ_vfc = [0.33, 0.33, 0.33, 0.33, 0.33]  # field capacity
    θ_vwp = [0.13, 0.13, 0.13, 0.13, 0.13]  # wilt point
    κ_dry = [6.3, 6.3, 6.3, 6.3, 6.3]  # thermal conductivity
    ψ_sat = [0.21, 0.22, 0.23, 0.24, 0.26]  # water potential at sat

  elseif SoilType == 6 # sandy clay loam
    b = [4.0, 4.2, 4.4, 4.6, 4.8]
    Ksat = [0.0000012, 0.00000108, 0.0000096, 0.0000072, 0.0000024]
    porosity = fill(0.398, 5)
    θ_vfc = fill(0.26, 5)
    θ_vwp = fill(0.15, 5)
    κ_dry = fill(7.0, 5)
    ψ_sat = [0.28, 0.29, 0.30, 0.31, 0.33]

  elseif SoilType == 7 # clay loam
    b = [5.2, 5.4, 5.6, 5.8, 6.0]
    Ksat = [0.00000064, 0.00000058, 0.00000051, 0.00000038, 0.00000013]
    porosity = fill(0.464, 5)
    θ_vfc = fill(0.32, 5)
    θ_vwp = fill(0.20, 5)
    κ_dry = [5.8, 5.8, 5.7, 5.8, 5.8]
    ψ_sat = [0.26, 0.27, 0.28, 0.29, 0.31]

  elseif SoilType == 8 # silty clay loam
    b = [6.6, 6.8, 7.0, 7.2, 7.4]
    Ksat = [0.00000042, 0.00000038, 0.00000034, 0.000000252, 0.000000084]
    porosity = fill(0.471, 5)
    θ_vfc = fill(0.37, 5)
    θ_vwp = fill(0.32, 5)
    κ_dry = fill(4.2, 5)
    ψ_sat = [0.33, 0.34, 0.35, 0.36, 0.38]

  elseif SoilType == 9 # sandy clay
    b = [6.0, 6.2, 6.4, 6.6, 6.8]
    Ksat = [0.00000033, 0.0000003, 0.000000264, 0.000000198, 0.000000066]
    porosity = fill(0.430, 5)
    θ_vfc = fill(0.34, 5)
    θ_vwp = fill(0.24, 5)
    κ_dry = fill(6.3, 5)
    ψ_sat = [0.29, 0.30, 0.31, 0.32, 0.34]

  elseif SoilType == 10 # silty clay
    b = [7.9, 8.1, 8.3, 8.5, 8.7]
    Ksat = [0.00000025, 0.000000225, 0.0000002, 0.00000015, 0.00000005]
    porosity = fill(0.479, 5)
    θ_vfc = fill(0.39, 5)
    θ_vwp = fill(0.25, 5)
    κ_dry = fill(4.0, 5)
    ψ_sat = [0.34, 0.35, 0.36, 0.37, 0.39]

  elseif SoilType == 11 # clay
    b = [7.6, 7.8, 8.0, 8.2, 8.4]
    Ksat = [0.00000017, 0.000000153, 0.000000136, 0.000000102, 0.000000034]
    porosity = fill(0.475, 5)
    θ_vfc = fill(0.40, 5)
    θ_vwp = fill(0.27, 5)
    κ_dry = fill(4.4, 5)
    ψ_sat = [0.37, 0.38, 0.39, 0.40, 0.42]

  else # default
    b = [7.6, 7.8, 8.0, 8.2, 8.4]
    Ksat = [0.00000017, 0.000000153, 0.000000136, 0.000000102, 0.000000034]
    porosity = fill(0.475, 5)
    θ_vfc = fill(0.40, 5)
    θ_vwp = fill(0.27, 5)
    κ_dry = fill(4.4, 5)
    ψ_sat = [0.37, 0.38, 0.39, 0.40, 0.42]
  end

  p.b[1:5] .= b
  p.Ksat[1:5] .= Ksat
  p.θ_sat[1:5] .= porosity
  p.θ_vfc[1:5] .= θ_vfc
  p.θ_vwp[1:5] .= θ_vwp
  p.κ_dry[1:5] .= κ_dry
  p.ψ_sat[1:5] .= ψ_sat
  return p
end


"""
update `f_root`
"""
function SoilRootFraction(soil::Soil)
  cum_depth = zeros(soil.n_layer)

  # For the first layer
  cum_depth[1] = soil.dz[1]
  soil.f_root[1] = 1 - (soil.r_root_decay)^(cum_depth[1] * 100)

  # For 2 to n_layer - 1
  for i in 2:(soil.n_layer-1)
    cum_depth[i] = cum_depth[i-1] + soil.dz[i]
    soil.f_root[i] = (soil.r_root_decay)^(cum_depth[i-1] * 100) - (soil.r_root_decay)^(cum_depth[i] * 100)
  end

  # For the last layer
  soil.f_root[soil.n_layer] = (soil.r_root_decay)^(cum_depth[soil.n_layer-1] * 100)
end


function Init_Soil_Status(p::Soil, Tsoil::Float64, Tair::Float64, θ0::Float64, snowdepth::Float64)
  d_t = clamp(Tsoil - Tair, -5.0, 5.0)

  # p.z_water = 0.0
  # p.r_rain_g = 0.0
  p.z_snow = snowdepth

  temp_scale_factors = [0.4, 0.5, 1.0, 1.2, 1.4]
  moisture_scale_factors = [0.8, 1.0, 1.05, 1.10, 1.15]

  for i in 1:p.n_layer
    p.Tsoil_c[i] = Tair + temp_scale_factors[i] * d_t
    p.Tsoil_p[i] = Tair + temp_scale_factors[i] * d_t
    p.θ[i] = moisture_scale_factors[i] * θ0
    p.θ_prev[i] = moisture_scale_factors[i] * θ0
  end

  for i in 1:p.n_layer
    if p.Tsoil_c[i] < -1.0
      p.ice_ratio[i] = 1.0
    elseif p.Tsoil_c[i] > 0
      p.ice_ratio[i] = 0
    else
      p.ice_ratio[i] = -p.Tsoil_c[i] / 1.0
    end
  end
end
