function besp_main(d::DataFrame, lai::Vector;
  lon::FT=120.0, lat::FT=20.0,
  VegType::Int=25, SoilType::Int=8,
  clumping::FT=0.85,
  Tsoil0::FT=2.2, θ0::FT=0.4115, z_snow0::FT=0.0,
  r_drainage::FT=0.5, r_root_decay::FT=0.95,
  version="julia", debug=false, fix_snowpack=true, kw...) where {FT<:AbstractFloat}

  meteo = Met()
  mid_res = Results()
  mid_ET = OutputET()
  Ra = Radiation()
  var = InterTempVars()

  theta = readVegParam(VegType)  # n = 48
  vegpar = theta2par(theta)
  theta = par2theta(vegpar; clumping, VegType) # 为移除readVegParam铺垫

  # coef = readcoef(VegType, SoilType) # n = 48, soil respiration module

  n = size(d, 1)
  vars = fieldnames(Results) |> collect
  vars_ET = fieldnames(OutputET) |> collect

  df_out = DataFrame(zeros(n, length(vars)), vars)
  df_ET = DataFrame(zeros(n, length(vars_ET)), vars_ET)

  # TODO: pass depth of soil as a parameter
  Tsoil = zeros(n, layer)
  θ = zeros(n, layer)

  if version == "julia"
    soil = Soil()
    state = State()
    state_n = State()
  elseif version == "c"
    soil = Soil_c()
    state = zeros(41)
    state_n = zeros(41)
  end
  Ta = d.tem[1]

  # (; r_drainage, r_root_decay) = vegpar
  Init_Soil_var_o(soil, state, Ta; VegType, SoilType,
    r_drainage, r_root_decay,
    Tsoil0, θ0, z_snow0
  )

  for i = 1:n
    jday = d.day[i]
    hour = d.hour[i]
    _day = ceil(Int, i / 24)
    mod(_day, 50) == 0 && (hour == 1) && println("Day = $_day")

    _lai = lai[_day] * theta[3] / clumping # re-calculate LAI & renew clump index
    fill_meteo!(meteo, d, i) # 驱动数据

    CosZs = s_coszs(jday, hour, lat, lon) # cos_solar zenith angle

    # /***** start simulation modules *****/
    if version == "julia"
      inter_prg_jl(jday, hour, _lai, clumping, vegpar, meteo, CosZs,
        state, soil,
        Ra, mid_res, mid_ET, var; VegType, debug, fix_snowpack)
    elseif version == "c"
      inter_prg_c(jday, hour, _lai, clumping, theta, meteo, CosZs,
        state, state_n, soil,
        Ra, mid_res, mid_ET, var; debug)
      state .= state_n # state variables
    end

    Tsoil[i, :] .= soil.Tsoil_c[1:layer]
    θ[i, :] .= soil.θ[1:layer]

    fill_res!(df_out, mid_res, i)
    fill_res!(df_ET, mid_ET, i)
  end

  df_out, df_ET, Tsoil, θ
end
