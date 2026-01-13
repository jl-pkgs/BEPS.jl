function beps_modern(d::DataFrame, lai::Vector;
  model::Union{Nothing,BEPSmodel}=nothing,
  lon::FT=120.0, lat::FT=20.0,
  VegType::Int=25, SoilType::Int=8,
  clumping::FT=0.85,
  Tsoil0::FT=2.2, θ0::FT=0.4115, z_snow0::FT=0.0,
  r_drainage::FT=0.5, r_root_decay::FT=0.95,
  debug=false, fix_snowpack=true, kw...) where {FT<:AbstractFloat}

  meteo = Met()
  mid_res = Results()
  mid_ET = OutputET()
  Ra = Radiation()
  var = TransientCache()

  if isnothing(model)
    theta_raw = readVegParam(VegType)  # n = 48
    vegpar = theta2par(theta_raw)
  else
    vegpar = model.veg
    (; r_drainage, r_root_decay) = model
  end

  n = size(d, 1)
  vars = fieldnames(Results) |> collect
  vars_ET = fieldnames(OutputET) |> collect

  df_out = DataFrame(zeros(n, length(vars)), vars)
  df_ET = DataFrame(zeros(n, length(vars_ET)), vars_ET)

  Tsoil = zeros(n, layer) # layer is defined in Constant.jl
  θ = zeros(n, layer)

  soil = Soil()
  state = State()
  Ta = d.tem[1]

  ## 土壤水力参数也应从model中读取
  Init_Soil_var(soil, state, Ta; VegType, SoilType,
    r_drainage=r_drainage, r_root_decay=r_root_decay,
    Tsoil0, θ0, z_snow0
  )
  !isnothing(model) && init_soil!(soil, model)

  for i = 1:n
    jday = d.day[i]
    hour = d.hour[i]

    _day = ceil(Int, i / 24)
    fill_meteo!(meteo, d, i)

    _lai = lai[_day]
    CosZs = s_coszs(jday, hour, lat, lon)

    # Main simulation step
    inter_prg_jl(jday, hour, _lai, clumping, vegpar, meteo, CosZs,
      state, soil,
      Ra, mid_res, mid_ET, var; VegType, debug, fix_snowpack)

    Tsoil[i, :] .= soil.Tsoil_c[1:layer]
    θ[i, :] .= soil.θ[1:layer]

    fill_res!(df_out, mid_res, i)
    fill_res!(df_ET, mid_ET, i)
  end
  df_out, df_ET, Tsoil, θ
end
