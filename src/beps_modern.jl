using ModelParams
using Statistics
using ModelParams: of_RMSE

function init_model()
  theta = ReadParamVeg(VegType)  # n = 48
  vegpar = theta2par(theta)
  theta = par2theta(vegpar; clumping, VegType) # 为移除ReadParamVeg铺垫
end

function besp_modern(d::DataFrame, lai::Vector; model::Union{Nothing,BEPSmodel}=nothing,
  lon::FT=120.0, lat::FT=20.0,
  VegType::Int=25, SoilType::Int=8, clumping::FT=0.85,
  Tsoil0::FT=2.2, θ0::FT=0.4115, z_snow0::FT=0.0,
  r_drainage::FT=0.5, r_root_decay::FT=0.95,
  debug=false, fix_snowpack=true, kw...) where {FT<:AbstractFloat}

  meteo = Met()
  mid_res = Results()
  mid_ET = OutputET()
  Ra = Radiation()
  var = TransientCache()

  if isnothing(model)
    theta = ReadParamVeg(VegType)  # n = 48
    vegpar = theta2par(theta)
    theta = par2theta(vegpar; clumping, VegType) # 为移除ReadParamVeg铺垫
  else
    vegpar = model.veg
    theta = par2theta(vegpar; clumping, VegType)
    (; r_drainage, r_root_decay) = model
  end

  n = size(d, 1)
  vars = fieldnames(Results) |> collect
  vars_ET = fieldnames(OutputET) |> collect

  ## OUTPUTs
  df_out = DataFrame(zeros(n, length(vars)), vars)
  df_ET = DataFrame(zeros(n, length(vars_ET)), vars_ET)

  output_Tsoil = zeros(n, layer) ## 返回变量
  output_θ = zeros(n, layer)

  soil = Soil()
  state = State()

  Ta = d.tem[1] # 第一天的温度
  Init_Soil_var(soil, state, Ta;
    VegType, SoilType, r_drainage, r_root_decay,
    Tsoil0, θ0, z_snow0
  )
  !isnothing(model) && Sync_Param_to_Soil!(soil, model) # 这是何意?

  for i = 1:n
    jday = d.day[i]
    hour = d.hour[i]
    _day = ceil(Int, i / 24)
    mod(_day, 50) == 0 && (hour == 1) && println("Day = $_day")

    _lai = lai[_day] * theta[3] / clumping # re-calculate LAI & renew clump index
    fill_meteo!(meteo, d, i) # 驱动数据

    CosZs = s_coszs(jday, hour, lat, lon) # cos_solar zenith angle

    inter_prg_jl(jday, hour, _lai, clumping, vegpar, meteo, CosZs,
      state, soil,
      Ra, mid_res, mid_ET, var; VegType, debug, fix_snowpack)

    output_Tsoil[i, :] .= soil.Tsoil_c[1:layer]
    output_θ[i, :] .= soil.θ[1:layer]
    fill_res!(df_out, mid_res, i)
    fill_res!(df_ET, mid_ET, i)
  end

  df_out, df_ET, output_Tsoil, output_θ
end
