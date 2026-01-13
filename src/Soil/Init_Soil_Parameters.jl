"""
    Init_Soil_Parameters(p::Soil, VegType::Integer, SoilType::Integer, r_root_decay::Float64)

Initialize soil parameters

- `Ksat`         : saturated hydraulic conductivity
- `porosity`     : porosity
- `θ_vfc`        : field capacity
- `θ_vwp`        : wilt point
- `ψ_sat`        : water potential at saturate
- `κ_dry`        : thermal conductivity
"""
function Init_Soil_Parameters(p::Soil, VegType::Integer, SoilType::Integer, r_root_decay::Float64)
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

  idx = (1 <= SoilType <= 11) ? SoilType : 11
  par = SOIL_PARAMS[idx]

  p.ρ_soil[1:5] .= [1300.0, 1500.0, 1517.0, 1517.0, 1517.0] # from flux tower.
  p.V_SOM[1:5] .= [5, 2, 1, 1, 0.3]

  p.b[1:5] .= par.b
  p.Ksat[1:5] .= par.K_sat
  p.θ_sat[1:5] .= fill(par.θ_sat, 5)
  p.θ_vfc[1:5] .= fill(par.θ_vfc, 5)
  p.θ_vwp[1:5] .= fill(par.θ_vwp, 5)
  p.κ_dry[1:5] .= fill(par.κ_dry, 5)
  p.ψ_sat[1:5] .= par.ψ_sat
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
