export InitState0, predict, loss, optim, BEPS_GOF

using ModelParams: sceua, of_KGE, of_NSE, GOF


function InitState0(model::ParamBEPS{FT}, forcing::MetSeries{FT}) where {FT<:AbstractFloat}
  Ta = forcing.Tair[1]
  Tsoil0 = Ta
  θ0 = model.hydraulic.θ_vfc[1]  # 初始化为田间持水量（避免低于凋萎含水量）
  z_snow0 = 0.0
  BEPS._init_state(model, Tsoil0, Ta, θ0, z_snow0) # state
end


function predict(theta::Vector{FT}, model::ParamBEPS{FT}, 
  forcing::MetSeries{FT}, lai::Vector{FT}, dates::Vector{DateTime}; 
  paths, lon::FT, lat::FT) where {FT<:AbstractFloat}

  model = deepcopy(model)
  params = parameters(model; paths)
  theta_prev = params.value
  BEPS.update!(model, paths, theta)

  state = InitState0(model, forcing)
  df_fluxes, df_ET, states, caches = simulate(forcing, lai, dates;
    ps=model, state, lon, lat)

  BEPS.update!(model, paths, theta_prev) # 恢复原参数值
  df_fluxes, df_ET, states, caches
end

# dates: UTC时间，需转换为当地时间
function loss(theta::Vector{FT}, model::ParamBEPS{FT}, 
  forcing::MetSeries{FT}, lai::Vector{FT}, dates::Vector{DateTime}; 
  paths, lon::FT, lat::FT, depths_SM::Vector{FT}, depths_TS::Vector{FT}, FluxDay::DataFrame) where {FT<:AbstractFloat}

  df_fluxes, df_ET, states, caches = predict(theta, model, forcing, lai, dates; 
    paths, lon, lat)

  dates_local = dates .+ Hour(8) # [UTC] -> [local]
  gof, data_sim, data_obs = BEPS_GOF(df_fluxes, states, dates_local, FluxDay;
    depths_SM, depths_TS)
  -mean(gof.Flux[1:2, :NSE])
end


# `dates`: 注意是UTC时间
"""
- `kw`: 
"""
function optim(model, forcing, lai, dates; paths, maxn=200, kw_loss...)
  _loss(theta) = loss(theta, model, forcing, lai, dates; paths, kw_loss...)

  params = parameters(model; paths)
  lb = map(x -> FT(x[1]), params.bound)
  ub = map(x -> FT(x[2]), params.bound)

  u0 = params.value
  theta, feval, exitflag = sceua(_loss, u0, lb, ub; maxn, verbose=true, parallel=true)
  theta
end



# GLOBAL: depths_SM, depths_TS
function BEPS_GOF(df_fluxes, states, dates_hour, FluxALL; depths_SM, depths_TS)
  ## 模拟
  (; lai, GPP_obs, ET_obs, Hs_obs) = FluxALL
  vars_TS = map(i -> Symbol("TS_$(Int(depths_TS[i] * 100))cm"), eachindex(depths_TS))
  vars_SM = map(i -> Symbol("SM_$(Int(depths_SM[i] * 100))cm"), eachindex(depths_SM))
  TS_obs = FluxALL[:, vars_TS] |> Matrix
  SM_obs = FluxALL[:, vars_SM] |> Matrix

  ## 观测
  (dates, GPP, ET, Hs) = agg_daily(df_fluxes, dates_hour)

  ## SM GOF
  SM_sim = states.vectors.θ
  TS_sim = states.vectors.Tsoil_c

  SM_day = agg_daily(SM_sim, dates_hour)
  TS_day = agg_daily(TS_sim, dates_hour)

  ## SM & TS
  SM = interp_depths(SM_day, depths_SM)
  TS = interp_depths(TS_day, depths_TS)

  gof_SM = map(i -> (; var=vars_SM[i], GOF(SM_obs[:, i], SM[:, i])...),
    eachindex(depths_SM)) |> DataFrame
  gof_TS = map(i -> (; var=vars_TS[i], GOF(TS_obs[:, i], TS[:, i])...),
    eachindex(depths_TS)) |> DataFrame

  gof_Flux = DataFrame([
    (; var="GPP", GOF(GPP_obs, GPP)...),
    (; var="ET", GOF(ET_obs, ET)...),
    (; var="Hs", GOF(Hs_obs, Hs)...)
  ])

  gof = (; Flux=gof_Flux, SM=gof_SM, TS=gof_TS)

  data_sim = hcat(
    DataFrame(; date=dates, GPP, ET, Hs),
    DataFrame(SM, vars_SM),
    DataFrame(TS, vars_TS))
  data_obs = hcat(
    DataFrame(; date=dates, GPP=GPP_obs, ET=ET_obs, Hs=Hs_obs),
    DataFrame(SM_obs, vars_SM),
    DataFrame(TS_obs, vars_TS))
  gof, data_sim, data_obs
end
