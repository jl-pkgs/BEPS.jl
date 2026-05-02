# see `VARS_SCALAR` and `VARS_VECTOR` for details
const DEFAULT_VARS_EXPORT = [
    :z_water, :ρ_snow, :z_snow, :r_rain_g, :f_soilwater,
    :θ, :Tsoil_c, :ETi]

function besp_main(forcing::MetSeries, lai::Vector, dates;
    lon::FT=120.0, lat::FT=20.0,
    VegType::Int=25, SoilType::Int=8, clumping::FT=0.85,
    Tsoil0::FT=2.2, θ0::FT=0.4115, z_snow0::FT=0.0,
    r_drainage::FT=0.5, r_root_decay::FT=0.95,
    VARS_EXPORT::Vector{Symbol}=DEFAULT_VARS_EXPORT,
    version="julia",
    fix_snowpack=true,
    fix_Ta_annual=true,
    kw...) where {FT<:AbstractFloat}

    ntime = forcing.ntime
    met = Met()
    mid_flux = Flux()
    mid_ET = ETFlux()
    Ra = Radiation()
    cache = LeafCache()

    fluxes = FluxSeries(; ntime)
    fluxes_ET = ETSeries(; ntime)
    SF, VF = split_vars(VARS_EXPORT)
    states = StateSeries(SF, VF, layer, ntime)
    # states = nothing

    ## 初始化参数和状态变量
    Ta = forcing.Tair[1]

    # 使用统一的 setup 函数初始化
    soil, state, params = setup_model(VegType, SoilType;
        version, Ta, Tsoil=Tsoil0, θ0, z_snow=z_snow0, r_drainage, r_root_decay, FT)
    state_n = deepcopy(state)

    theta = par2theta(params.veg; clumping, VegType)

    Ta_annual = mean(forcing.Tair)
    jdays = dayofyear.(dates)
    hours = hour.(dates) .+ 1  # 转为1-based (1:24)

    for i = 1:ntime
        # add a progress
        jday = jdays[i]
        hour = hours[i]

        fill_met!(met, forcing, i) # 驱动数据
        k = ceil(Int, i / 24)
        _lai = lai[k] * theta[3] / clumping # re-calculate LAI & renew clump index

        # /***** start simulation modules *****/
        if version == "julia"
            inter_prg_jl(jday, hour, lon, lat, _lai, clumping,
                Ra, met, params, state, mid_flux, mid_ET, cache; fix_Ta_annual, fix_snowpack, Ta_annual)
            save_state!(states, state, i, SF, VF)
        elseif version == "c"
            inter_prg_c(jday, hour, lon, lat, _lai, clumping,
                Ra, met, theta, state, state_n, soil, mid_flux, mid_ET, cache;)
            state .= state_n # state variables
        end

        fluxes[i] = mid_flux
        fluxes_ET[i] = mid_ET
    end
    DataFrame(fluxes), DataFrame(fluxes_ET), states
end
