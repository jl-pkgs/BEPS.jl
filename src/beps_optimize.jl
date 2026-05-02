function beps_optimize(forcing::MetSeries, lai::Vector, dates::AbstractVector,
  model::ParamBEPS, obs::AbstractVector;
  col_sim::Symbol=:ET,
  # SCE-UA 优化器参数
  maxn=2000, kstop=5, f_reltol=0.0001, x_reltol=0.0001, seed=1,
  n_complex=5, size_complex=nothing, size_simplex=nothing, n_evolu=nothing,
  n_pop=nothing, verbose=false, parallel=true,
  kwargs...)

  # state0 初始条件与优化参数无关，只初始化一次；beps_modern 内部会 deepcopy
  Ta = Float64(forcing.Tair[1])
  state0, _ = setup(model; Ta)
  x0, lb, ub, paths = get_opt_info(model)

  function cal_func(x)
    m = deepcopy(model)
    update!(m, paths, x)
    try
      df_out, df_ET, _ = beps_modern(forcing, lai, dates; ps=m, state=state0, kwargs...)
      sim = col_sim ∈ propertynames(df_out) ? df_out[!, col_sim] : df_ET[!, col_sim]
      return of_RMSE(obs, sim)
    catch e
      e isa DomainError || @warn "beps_optimize: unexpected error" exception=e
      return Inf  # 参数违反物理约束时返回大值
    end
  end

  sceua_kw = (; maxn, kstop, f_reltol, x_reltol, seed, n_complex, verbose, parallel)
  opt_extras = (; size_complex, size_simplex, n_evolu, n_pop)
  sceua_kw = merge(sceua_kw,
    NamedTuple(k => v for (k, v) in pairs(opt_extras) if !isnothing(v)))
  bestx, bestf, exitflag = sceua(cal_func, x0, lb, ub; sceua_kw...)

  update!(model, paths, bestx)
  return model, bestf
end
