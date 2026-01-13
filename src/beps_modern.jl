function besp_modern(d::DataFrame, lai::Vector; model::Union{Nothing,BEPSmodel}=nothing,
  lon::FT=120.0, lat::FT=20.0,
  VegType::Int=25, clumping::FT=0.85,
  Tsoil0::FT=2.2, θ0::FT=0.4115, z_snow0::FT=0.0,
  debug=false, fix_snowpack=true, kw...) where {FT<:AbstractFloat}

  meteo = Met()
  mid_res = Results()
  mid_ET = OutputET()
  Ra = Radiation()
  cache = TransientCache()

  ntime = size(d, 1)
  vars = fieldnames(Results) |> collect
  vars_ET = fieldnames(OutputET) |> collect

  ## OUTPUTs
  df_out = DataFrame(zeros(ntime, length(vars)), vars)
  df_ET = DataFrame(zeros(ntime, length(vars_ET)), vars_ET)
  output_Tsoil = zeros(ntime, layer) ## 返回变量
  output_θ = zeros(ntime, layer)

  Ta = d.tem[1] # 第时刻的温度
  state = State()
  soil = Soil()
  Sync_Param_to_Soil!(soil, model)
  Init_Soil_T_θ!(soil, Tsoil0, Ta, θ0, z_snow0)
  Sync_Soil_to_State!(soil, state, Ta)

  vegpar = model.veg

  for i = 1:ntime
    jday = d.day[i]
    hour = d.hour[i]
    _day = ceil(Int, i / 24)
    mod(_day, 50) == 0 && (hour == 1) && println("Day = $_day")

    _lai = lai[_day]
    fill_meteo!(meteo, d, i) # 驱动数据

    CosZs = s_coszs(jday, hour, lat, lon) # cos_solar zenith angle

    inter_prg_jl(jday, hour, _lai, clumping, vegpar, meteo, CosZs,
      state, soil,
      Ra, mid_res, mid_ET, cache; VegType, debug, fix_snowpack)

    output_Tsoil[i, :] .= soil.Tsoil_c[1:layer]
    output_θ[i, :] .= soil.θ[1:layer]
    fill_res!(df_out, mid_res, i)
    fill_res!(df_ET, mid_ET, i)
  end

  df_out, df_ET, output_Tsoil, output_θ
end
