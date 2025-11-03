@bounds @with_kw mutable struct BEPSmodel2{FT<:AbstractFloat}
  N::Int = 5
  r_drainage::FT = Cdouble(0.50) | (0.2, 0.7)      # ? 地表排水速率（地表汇流），可考虑采用曼宁公式
  r_root_decay::FT = Cdouble(0.95) | (0.85, 0.999) # ? 根系分布衰减率, decay_rate_of_root_distribution

  ψ_min::FT = Cdouble(33.0)  # * 气孔关闭对应水势，33kPa，可根据植被类型指定
  alpha::FT = Cdouble(0.4)   # * 土壤水限制因子参数，He 2017 JGR-B, Eq. 4

  hydraulic::ParamSoilHydraulicLayers{FT} = ParamSoilHydraulicLayers{FT,N}()
  thermal::ParamSoilThermalLayers{FT} = ParamSoilThermalLayers{FT,N}()

  veg::VegParam{FT} = VegParam{FT}()
end


function BEPSmodel(; VegType::Int=25, SoilType::Int=8,
  r_drainage::FT=0.5, r_root_decay::FT=0.95,
) where {FT<:AbstractFloat}
  theta = readVegParam(VegType)  # n = 48
  vegpar = theta2par(theta)

  soil = Soil()
  state = State()

  # (; r_drainage, r_root_decay) = vegpar
  # 初始化state, 同时填充params
  Init_Soil_var_o(soil, state, Ta; VegType, SoilType,
    r_drainage, r_root_decay,
    Tsoil0, θ0, z_snow0
  )
end


function beps_modern(d::DataFrame, lai::Vector;
  lon::FT=120.0, lat::FT=20.0,
  VegType::Int=25, SoilType::Int=8,
  # clumping::FT=0.85,
  # Tsoil0::FT=2.2, θ0::FT=0.4115, z_snow0::FT=0.0,
  debug=false, fix_snowpack=true, kw...) where {FT<:AbstractFloat}

  meteo = Met()
  mid_res = Results()
  mid_ET = OutputET()
  Ra = Radiation()
  var = InterTempVars()

  n = size(d, 1)
  vars = fieldnames(Results) |> collect
  vars_ET = fieldnames(OutputET) |> collect

  df_out = DataFrame(zeros(n, length(vars)), vars)
  df_ET = DataFrame(zeros(n, length(vars_ET)), vars_ET)

  Tsoil = zeros(n, layer) # soil temperature
  θ = zeros(n, layer)     # soil moisture
  Ta = d.tem[1]

  for i = 1:n
    jday = d.day[i]
    hour = d.hour[i]

    _day = ceil(Int, i / 24) # progress
    mod(_day, 50) == 0 && (hour == 1) && println("Day = $_day")

    # LAI是按天输入的
    _lai = lai[_day] # * theta[3] / clumping # re-calculate LAI & renew clump index, 没道理
    fill_meteo!(meteo, d, i) # 驱动数据

    CosZs = s_coszs(jday, hour, lat, lon) # cos_solar zenith angle, UTC time

    # /***** start simulation modules *****/
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
