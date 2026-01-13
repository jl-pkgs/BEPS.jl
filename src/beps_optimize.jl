function beps_optimize(d::DataFrame, lai::Vector, model::BEPSmodel, obs::AbstractVector;
  col_sim::Symbol=:ET,
  maxn=2000,
  kstop=5, pcento=0.01, peps=0.001, iseed=1, iniflg=0,
  kwargs...)

  x0, lb, ub, paths = get_opt_info(model)

  function cal_func(x)
    update!(model, paths, x)
    df_out, df_ET, _, _ = besp_main(d, lai; model=model, kwargs...)

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
    maxn=maxn, kstop=kstop, pcento=pcento, peps=peps, iseed=iseed, iniflg=iniflg)

  update!(model, paths, bestx)
  return model, bestf
end
