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
function Init_Soil_Parameters(soil::Soil, VegType::Integer, SoilType::Integer, r_root_decay::Float64)
  soil.n_layer = 5

  if VegType == 6 || VegType == 9 # DBF or EBF, low constaint threshold
    soil.ψ_min = 10.0 # ψ_min
    soil.alpha = 1.5
  else
    soil.ψ_min = 33.0 # ψ_min
    soil.alpha = 0.4
  end

  soil.dz[1:5] .= [0.05, 0.10, 0.20, 0.40, 1.25] # BEPS V2023
  # z = [0, 5, 15, 25, 35, 45, 55.0] ./ 100
  # z_mid = (z[1:end-1] .+ z[2:end]) ./ 2
  # dz = diff(z)
  # n = length(z)
  # p.dz[1:n-1] = dz
  soil.r_root_decay = r_root_decay
  UpdateRootFraction!(soil)

  idx = (1 <= SoilType <= 11) ? SoilType : 11
  par = SOIL_PARAMS[idx]

  soil.ρ_soil[1:5] .= [1300.0, 1500.0, 1517.0, 1517.0, 1517.0] # from flux tower.
  soil.V_SOM[1:5] .= [5, 2, 1, 1, 0.3] # volume fraction, 0-1; bug 20260113, unit error

  soil.b[1:5] .= par.b
  soil.Ksat[1:5] .= par.K_sat
  soil.θ_sat[1:5] .= fill(par.θ_sat, 5)
  soil.θ_vfc[1:5] .= fill(par.θ_vfc, 5)
  soil.θ_vwp[1:5] .= fill(par.θ_vwp, 5)
  soil.κ_dry[1:5] .= fill(par.κ_dry, 5)
  soil.ψ_sat[1:5] .= par.ψ_sat
  return soil
end


# landcover::Int, soil_type::Int, Tsoil, soilwater, snowdepth
# initialize soil conditions, read soil parameters and set depth
function Init_Soil_var(soil::AbstractSoil, state::Union{State,Vector}, Ta::FT;
  VegType::Int=25, SoilType::Int=8,
  r_drainage::FT, r_root_decay::FT,
  Tsoil0::FT=Ta, θ0::FT=0.2, z_snow0::FT=0.0,
  ignored...
) where {FT<:Real}
  # r_drainage = param[27]
  # r_root_decay = param[28]
  soil.r_drainage = r_drainage
  Init_Soil_Parameters(soil, VegType, SoilType, r_root_decay)
  Init_Soil_T_θ!(soil, Tsoil0, Ta, θ0, z_snow0)
  Sync_Soil_to_State!(soil, state, Ta)
end


export Init_Soil_var
