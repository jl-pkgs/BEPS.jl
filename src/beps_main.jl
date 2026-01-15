function besp_main(d::DataFrame, lai::Vector; params::Union{Nothing,ParamBEPS}=nothing,
  lon::FT=120.0, lat::FT=20.0,
  VegType::Int=25, SoilType::Int=8, clumping::FT=0.85,
  Tsoil0::FT=2.2, θ0::FT=0.4115, z_snow0::FT=0.0,
  r_drainage::FT=0.5, r_root_decay::FT=0.95,
  version="julia", fix_snowpack=true, kw...) where {FT<:AbstractFloat}

  met = Met()
  mid_res = Results()
  mid_ET = OutputET()
  Ra = Radiation()
  cache = LeafCache()

  n = size(d, 1)
  vars = fieldnames(Results) |> collect
  vars_ET = fieldnames(OutputET) |> collect

  df_out = DataFrame(zeros(n, length(vars)), vars)
  df_ET = DataFrame(zeros(n, length(vars_ET)), vars_ET)

  output_Tsoil = zeros(n, layer) ## 返回变量
  output_θ = zeros(n, layer)

  ## 初始化参数和状态变量
  Ta = d.tem[1]
  ps_veg = isnothing(params) ? InitParam_Veg(VegType; FT) : params.veg
  theta = par2theta(ps_veg; clumping, VegType)

  soil = version == "julia" ? Soil() : Soil_c()
  soil.r_drainage = r_drainage
  Init_Soil_Parameters(soil, VegType, SoilType, r_root_decay)
  Init_Soil_T_θ!(soil, Tsoil0, Ta, θ0, z_snow0)

  state = version == "julia" ? StateBEPS(soil) : zeros(41)
  state_n = deepcopy(state)
  InitState!(soil, state, Ta) # initialize state variables, for C version
  Params2Soil!(soil, params)  # put params into soil struct

  if isnothing(params)
    params = ParamBEPS(VegType, SoilType; FT=FT)
    version == "julia" && Soil2Params!(params, soil)
  end

  for i = 1:n
    jday = d.day[i]
    hour = d.hour[i]
    CosZs = s_coszs(jday, hour, lat, lon) # cos_solar zenith angle

    _day = ceil(Int, i / 24)
    mod(_day, 50) == 0 && (hour == 1) && println("Day = $_day")

    _lai = lai[_day] * theta[3] / clumping # re-calculate LAI & renew clump index
    fill_met!(met, d, i) # 驱动数据

    # /***** start simulation modules *****/
    if version == "julia"
      inter_prg_jl(jday, hour, CosZs, Ra, _lai, clumping, 
        met, params, state, mid_res, mid_ET, cache; fix_snowpack)
    elseif version == "c"
      inter_prg_c(jday, hour, CosZs, Ra, _lai, clumping, 
        met, theta, state, state_n, soil, mid_res, mid_ET, cache;)
      state .= state_n # state variables
    end

    output_Tsoil[i, :] .= soil.Tsoil_c[1:layer]
    output_θ[i, :] .= soil.θ[1:layer]

    fill_res!(df_out, mid_res, i)
    fill_res!(df_ET, mid_ET, i)
  end
  df_out, df_ET, output_Tsoil, output_θ
end
