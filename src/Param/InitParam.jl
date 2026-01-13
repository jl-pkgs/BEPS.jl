using JSON

const path_veg = joinpath(@__DIR__, "params", "ParamVeg.json")
const path_gen = joinpath(@__DIR__, "params", "ParamGeneral.json")

_types = ["ENF", "DNF", "DBF", "EBF", "Shrub-SH", "C4", "default"]
_codes = [1, 2, 6, 9, 13, 40, -1]

# const RTIMES = 24.0 # 呼吸作用系数转换常数 (day -> hour)

readGeneralParam() = JSON.parsefile(path_gen)

function readVegRaw(lc::Int=1)
  veg_data = JSON.parsefile(path_veg)
  type_idx = findfirst(x -> x == lc, _codes)
  type_str = type_idx !== nothing ? _types[type_idx] : "default"
  return veg_data[type_str]
end


"""
    InitParam_Veg(lc::Int=1; FT=Float64)

读取 JSON 配置文件并返回 VegParam 结构体。
"""
function InitParam_Veg(lc::Int=1; FT=Float64)
  veg_data = JSON.parsefile(path_veg)
  gen_data = JSON.parsefile(path_gen)

  type_idx = findfirst(x -> x == lc, _codes)
  type_str = type_idx !== nothing ? _types[type_idx] : "default"
  v = veg_data[type_str]

  return VegParam{FT}(
    LAI_max_o    = FT(v["LAI_max_o"]),
    LAI_max_u    = FT(v["LAI_max_u"]),
    α_canopy_vis = FT(v["albedo_canopy_vis"]),
    α_canopy_nir = FT(v["albedo_c_nir_observed"]),
    α_soil_sat   = FT(gen_data["albedo_saturated_soil"]),
    α_soil_dry   = FT(gen_data["albedo_dry_soil"]),
    z_canopy_o   = FT(v["z_canopy_o"]),
    z_canopy_u   = FT(v["z_canopy_u"]),
    z_wind       = FT(gen_data["the_height_to_measure_wind_speed"]),
    g1_w         = FT(v["g1_w"]),
    g0_w         = FT(gen_data["intercept_for_H2O_ball_berry"]),
    VCmax25      = FT(v["VCmax25"]),
    N_leaf       = FT(v["N_leaf"]),
    slope_Vc     = FT(v["slope_Vc"])
  )
end


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
function InitParam_Soil(p::Soil, VegType::Integer, SoilType::Integer, r_root_decay::Float64)
  p.n_layer = 5
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

# 这可以切分为: 

# if VegType == 6 || VegType == 9 # DBF or EBF, low constaint threshold
#   p.ψ_min = 10.0 # ψ_min
#   p.alpha = 1.5
# else
#   p.ψ_min = 33.0 # ψ_min
#   p.alpha = 0.4
# end
