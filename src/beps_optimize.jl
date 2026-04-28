function beps_optimize(d::DataFrame, lai::Vector, model::ParamBEPS, obs::AbstractVector;
  col_sim::Symbol=:ET,
  maxn=2000, kstop=5, f_reltol=0.0001, x_reltol=0.0001, seed=1,
  kwargs...)

  x0, lb, ub, paths = get_opt_info(model)

  function cal_func(x)
    # deepcopy 避免迭代间状态污染
    m = deepcopy(model)
    update!(m, paths, x)
    df_out = df_ET = nothing
    try
      df_out, df_ET, _, _ = beps_modern(d, lai; model=m, kwargs...)
    catch e
      e isa DomainError || @warn "beps_optimize: unexpected error" exception=e
      return Inf  # 参数违反物理约束时返回大值
    end

    sim = if col_sim in propertynames(df_out)
      df_out[!, col_sim]
    elseif col_sim in propertynames(df_ET)
      df_ET[!, col_sim]
    else
      error("Column $col_sim not found in output")
    end
    of_RMSE(obs, sim)
  end

  bestx, bestf, exitflag = sceua(cal_func, x0, lb, ub;
    maxn, kstop, f_reltol, x_reltol, seed)

  update!(model, paths, bestx)
  return model, bestf
end
