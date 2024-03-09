function besp_main(d::DataFrame, lai::Vector, par::NamedTuple;
  version="julia", debug=false, kw...)

  
  meteo = Met()
  mid_res = Results()
  mid_ET = OutputET()
  Ra = Radiation()
  var = InterTempVars()

  param = readparam(par.landcover)      # n = 48
  # coef = readcoef(par.landcover, par.soil_type) # n = 48, soil respiration module

  n = size(d, 1)
  vars = fieldnames(Results) |> collect
  vars_ET = fieldnames(OutputET) |> collect

  df_out = DataFrame(zeros(n, length(vars)), vars)
  df_ET = DataFrame(zeros(n, length(vars_ET)), vars_ET)

  # TODO: pass depth of soil as a parameter
  Tg = zeros(n, layer)
  θ = zeros(n, layer)

  var_o = zeros(41)
  var_n = zeros(41)
  v2last = zeros(41)

  if version == "julia"
    fun = inter_prg_jl
    p_soil = Soil()
  elseif version == "c"
    fun = inter_prg_c
    p_soil = Soil_c()
  end

  init = true
  
  for i = 1:n
    jday = d.day[i]
    hour = d.hour[i]

    _day = ceil(Int, i / 24)
    mod(_day, 50) == 0 && (hour == 1) && println("Day = $_day")

    _lai = lai[_day] * param[3] / par.clumping # re-calculate LAI & renew clump index

    fill_meteo!(meteo, d, i)

    if init
      Init_Soil_var_o(p_soil, var_o, meteo, param, par) # update p_soil and var_o
      init = false
    else
      var_o .= v2last
    end

    CosZs = s_coszs(jday, hour, par.lat, par.lon) # cos_solar zenith angle

    # /***** start simulation modules *****/
    fun(jday, hour, _lai, par.clumping, param, meteo, CosZs, var_o, var_n, p_soil,
      Ra, mid_res, mid_ET, var; debug)

    Tg[i, :] .= p_soil.Tsoil_c[1:layer]
    θ[i, :] .= p_soil.θ[1:layer]
    # Store updated variables array in temp array
    v2last .= var_n

    fill_res!(df_out, mid_res, i)
    fill_res!(df_ET, mid_ET, i)
  end

  df_out, df_ET, Tg, θ
end
