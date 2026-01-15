function besp_modern(d::DataFrame, lai::Vector; model::Union{Nothing,BEPSmodel}=nothing,
  lon::FT=120.0, lat::FT=20.0, clumping::FT=0.85,
  Tsoil0::FT=2.2, θ0::FT=0.4115, z_snow0::FT=0.0,
  fix_snowpack=true, kw...) where {FT<:AbstractFloat}

  met = Met()
  mid_res = Results()
  mid_ET = OutputET()
  Ra = Radiation()
  cache = LeafCache()

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
  Params2Soil!(soil, model)
  Init_Soil_T_θ!(soil, Tsoil0, Ta, θ0, z_snow0)
  InitState!(soil, state, Ta)

  ps_veg = model.veg

  for i = 1:ntime
    jday = d.day[i]
    hour = d.hour[i]
    CosZs = s_coszs(jday, hour, lat, lon) # cos_solar zenith angle
    
    _day = ceil(Int, i / 24)
    mod(_day, 50) == 0 && (hour == 1) && println("Day = $_day")

    _lai = lai[_day]
    fill_met!(met, d, i) # 驱动数据

    inter_prg_jl(jday, hour, CosZs, Ra, _lai, clumping, ps_veg, met,
      state, soil, mid_res, mid_ET, cache; fix_snowpack)

    output_Tsoil[i, :] .= soil.Tsoil_c[1:layer]
    output_θ[i, :] .= soil.θ[1:layer]
    fill_res!(df_out, mid_res, i)
    fill_res!(df_ET, mid_ET, i)
  end

  df_out, df_ET, output_Tsoil, output_θ
end
