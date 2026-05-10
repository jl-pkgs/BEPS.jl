using Plots
gr(; framestyle=:box)

function plot_result(sim, obs)
  dates = sim.date
  p1 = plot(dates, obs.GPP, label="GPP_obs")
  plot!(p1, dates, sim.GPP, label="GPP_sim")

  p2 = plot(dates, obs.ET, label="ET_obs")
  plot!(p2, dates, sim.ET, label="ET_sim")

  p3 = plot(dates, obs.Hs, label="Hs_obs")
  plot!(p3, dates, sim.Hs, label="Hs_sim")

  p = plot(p1, p2, p3, layout=(3, 1), size=(1400, 900))
  savefig(p, "Figure1_Fluxes.png")

  ## SM
  select_vars(df, starts="SM") = df[:, Cols(startswith(starts))] |> Matrix

  p = plot(dates, select_vars(obs, "SM"); label="OBS", size=(1400, 400))
  plot!(p, dates, select_vars(sim, "SM"); label="SIM")
  savefig(p, "Figure1_SM.png")

  p = plot(dates, select_vars(obs, "TS"); label="OBS", size=(1400, 400))
  plot!(p, dates, select_vars(sim, "TS"); label="SIM")
  savefig(p, "Figure1_TS.png")
end
