# DBF or EBF, low constaint threshold
function Sync_Param_to_Soil!(soil::Soil, model::BEPSmodel{FT}; BF=false) where {FT}
  soil.ψ_min = BF ? 10.0 : 33.0 # [m], about 0.10~0.33 MPa开始胁迫点
  soil.alpha = BF ? 1.5 : 0.4   # He 2017 JGR-B, Eq. 4

  (; hydraulic, thermal, N) = model
  soil.n_layer = Cint(N)
  soil.dz[1:5] .= [0.05, 0.10, 0.20, 0.40, 1.25] # BEPS V2023, 土壤层厚度[m]

  soil.r_drainage = Cdouble(model.r_drainage)
  soil.r_root_decay = Cdouble(model.veg.r_root_decay)
  UpdateRootFraction!(soil) # 更新根系分布

  soil.ψ_min = Cdouble(model.ψ_min)
  soil.alpha = Cdouble(model.alpha)

  soil.θ_vfc[1:N] .= Cdouble.(hydraulic.θ_vfc)
  soil.θ_vwp[1:N] .= Cdouble.(hydraulic.θ_vwp)
  soil.θ_sat[1:N] .= Cdouble.(hydraulic.θ_sat)
  soil.Ksat[1:N] .= Cdouble.(hydraulic.K_sat)
  soil.ψ_sat[1:N] .= Cdouble.(hydraulic.ψ_sat)
  soil.b[1:N] .= Cdouble.(hydraulic.b)

  soil.κ_dry[1:N] .= Cdouble.(thermal.κ_dry)
  soil.ρ_soil[1:N] .= Cdouble.(thermal.ρ_soil)
  soil.V_SOM[1:N] .= Cdouble.(thermal.V_SOM)
  return soil
end


# for (i=3;i<=8;i++)   var_o[i] = tem;
# for(i=9;i<=14;i++) var_o[i] = soil->Tsoil_p[i-9];
# for(i=21;i<=26;i++) var_o[i] = soil->θ_prev[i-21];
# for(i=27;i<=32;i++) var_o[i] = soil->ice_ratio[i-27];
function Sync_Soil_to_State!(soil::AbstractSoil, state::Vector, Ta)
  state .= 0
  for i = 1:6
    state[i+3] = Ta
    state[i+9] = soil.Tsoil_p[i]
    state[i+21] = soil.θ_prev[i]
    state[i+27] = soil.ice_ratio[i]
  end
  return nothing
end

function Sync_Soil_to_State!(soil::AbstractSoil, state::State, Ta)
  state.Ts .= Ta
  state.Ts_prev .= soil.Tsoil_p[1:5]
  state.θ_prev .= soil.θ_prev[1:5]
  state.ice_ratio .= soil.ice_ratio[1:5]
  return nothing
end
