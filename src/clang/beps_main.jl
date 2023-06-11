include("BEPS_helper.jl")


function besp_main(d::DataFrame, lai::Vector, par::NamedTuple; 
    version="julia", kw...)
  
  p_soil = Soil()
  meteo = ClimateData()
  mid_res = Results()
  mid_ET = OutputET()
  
  parameter = readparam(par.landcover)      # n = 48
  # coef = readcoef(par.landcover, par.soil_type) # n = 48, soil respiration module

  n = size(d, 1)
  vars = fieldnames(Results) |> collect
  vars_ET = fieldnames(OutputET) |> collect
  
  df_out = DataFrame(zeros(n, length(vars)), vars)
  df_ET = DataFrame(zeros(n, length(vars_ET)), vars_ET)
  
  var_o = zeros(41)
  var_n = zeros(41)
  v2last = zeros(41)

  # if version == "julia"
    fun = inter_prg_jl
  # elseif version == "c"
  #   fun = inter_prg
  # end

  for jday = 1:365
    if mod(jday, 50) == 0
      @show jday
    end
    _lai = lai[jday] * parameter[3] / par.clumping # re-calculate LAI & renew clump index

    for rstep = 1:24
      k = (jday - 1) * 24 + rstep
      flag = (jday == 1 && rstep == 1) ? 0 : 1

      fill_meteo!(meteo, d, k)

      if (flag == 0)
        init_soil!(p_soil, var_o, meteo, parameter, par) # update p_soil and var_o
      else
        var_o .= v2last
      end

      CosZs = s_coszs(jday, rstep - 1, par.lat, par.lon) # cos_solar zenith angle
      # /***** start simulation modules *****/
      
      fun(jday, rstep - 1, _lai, par.clumping, parameter, meteo, CosZs, var_o, var_n, p_soil, mid_res, mid_ET)
      # inter_prg_jl(jday, rstep - 1, _lai, par.clumping, parameter, meteo, CosZs, var_o, var_n, p_soil, mid_res)
      # Store updated variables array in temp array
      v2last .= var_n

      fill_res!(df_out, mid_res, k)
      fill_res!(df_ET, mid_ET, k)
    end # End of hourly loop
  end # End of daily loop
  df_out, df_ET
end
