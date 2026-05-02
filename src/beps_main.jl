function besp_main(d::DataFrame, lai::Vector;
  lon::FT=120.0, lat::FT=20.0,
  VegType::Int=25, SoilType::Int=8, clumping::FT=0.85,
  Tsoil0::FT=2.2, θ0::FT=0.4115, z_snow0::FT=0.0,
  r_drainage::FT=0.5, r_root_decay::FT=0.95,
  use_lrad::Bool=false,
  version="julia", fix_snowpack=true,
  verbose=true, kw...) where {FT<:AbstractFloat}

  d = standardize_forcing_columns(d)
  ntime = size(d, 1)

  met = Met()
  mid_flux = Flux()
  mid_ET = ETFlux()
  Ra = Radiation()
  cache = LeafCache()

  fluxes = FluxSeries(; ntime)
  fluxes_ET = ETSeries(; ntime)

  output_Tsoil = zeros(ntime, layer) ## 返回变量
  output_θ = zeros(ntime, layer)

  ## 初始化参数和状态变量
  Ta = d.Tair[1]

  # 使用统一的 setup 函数初始化
  soil, state, params = setup_model(VegType, SoilType;
    version, Ta, Tsoil=Tsoil0, θ0, z_snow=z_snow0, r_drainage, r_root_decay, FT)
  state_n = deepcopy(state)

  theta = par2theta(params.veg; clumping, VegType)

  for i = 1:ntime
    jday = d.day[i]
    hour = d.hour[i]
    CosZs = s_coszs(jday, hour, lat, lon) # cos_solar zenith angle

    _day = ceil(Int, i / 24)
    (mod(_day, 50) == 0 && (hour == 1) && verbose) && println("Day = $_day")

    _lai = lai[_day] * theta[3] / clumping # re-calculate LAI & renew clump index
    fill_met!(met, d, i; use_lrad) # 驱动数据

    # /***** start simulation modules *****/
    if version == "julia"
      inter_prg_jl(jday, hour, CosZs, Ra, _lai, clumping,
        met, params, state, mid_flux, mid_ET, cache; fix_snowpack)
    elseif version == "c"
      inter_prg_c(jday, hour, CosZs, Ra, _lai, clumping,
        met, theta, state, state_n, soil, mid_flux, mid_ET, cache;)
      state .= state_n # state variables
    end

    output_Tsoil[i, :] .= soil.Tsoil_c[1:layer]
    output_θ[i, :] .= soil.θ[1:layer]
    fluxes[i] = mid_flux
    fluxes_ET[i] = mid_ET
  end
  DataFrame(fluxes), DataFrame(fluxes_ET), output_Tsoil, output_θ
end
