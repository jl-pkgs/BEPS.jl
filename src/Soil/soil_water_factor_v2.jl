# Function to compute soil water stress factor
function soil_water_factor_v2(p::Soil)
  # @unpack ψ, ψ_sat, 
  (; ψ_min) = p
  θ = p.θ
  n = p.n_layer

  t1 = -0.02
  t2 = 2.0

  if p.ψ[1] <= 0.000001
    for i in 1:n
      p.ψ[i] = cal_ψ(θ[i], p.θ_sat[i], p.ψ_sat[i], p.b[i])
    end
  end

  for i in 1:n
    # psi_sr in m H2O! This is the old version. LHE; He 2017 JGR-B, Eq. 4
    p.fpsisr[i] = p.ψ[i] > ψ_min ? 1.0 / (1 + ((p.ψ[i] - ψ_min) / ψ_min)^p.alpha) : 1.0

    p.ft[i] = p.Tsoil_p[i] > 0.0 ? 1.0 - exp(t1 * p.Tsoil_p[i]^t2) : 0
    p.fpsisr[i] *= p.ft[i]
  end

  for i in 1:n
    p.dtt[i] = FW_VERSION == 1 ? p.f_root[i] * p.fpsisr[i] : p.f_root[i]
  end
  dtt_sum = sum(p.dtt) # 每层的土壤水分限制因子

  if dtt_sum < 0.000001
    p.f_soilwater = 0.1
  else
    fpsisr_sum = 0
    for i in 1:n
      p.dt[i] = p.dtt[i] / dtt_sum
      isnan(p.dt[i]) && println(p.dt[1])

      fpsisr_sum += p.fpsisr[i] * p.dt[i]
    end
    p.f_soilwater = max(0.1, fpsisr_sum)
  end
end
