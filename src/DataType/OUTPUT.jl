@with_kw mutable struct Results
  gpp_o_sunlit::Cdouble = 0.0
  gpp_u_sunlit::Cdouble = 0.0
  gpp_o_shaded::Cdouble = 0.0
  gpp_u_shaded::Cdouble = 0.0

  plant_resp::Cdouble = 0.0
  npp_o::Cdouble = 0.0
  npp_u::Cdouble = 0.0
  GPP::Cdouble = 0.0
  NPP::Cdouble = 0.0
  NEP::Cdouble = 0.0
  soil_resp::Cdouble = 0.0
  Net_Rad::Cdouble = 0.0

  SH::Cdouble = 0.0
  LH::Cdouble = 0.0
  Trans::Cdouble = 0.0
  Evap::Cdouble = 0.0

  z_water::Cdouble = 0.0
  z_snow::Cdouble = 0.0
  ρ_snow::Cdouble = 0.0
end


@with_kw mutable struct OutputET
  Trans_o::Cdouble = 0.0
  Trans_u::Cdouble = 0.0
  Eil_o::Cdouble = 0.0     # Ei of liquid
  Eil_u::Cdouble = 0.0 
  EiS_o::Cdouble = 0.0     # Ei of solid
  EiS_u::Cdouble = 0.0
  Evap_soil::Cdouble = 0.0 
  Evap_SW::Cdouble = 0.0   # evaporation from water pond
  Evap_SS::Cdouble = 0.0   # evaporation from snow pack
  Qhc_o::Cdouble = 0.0
  Qhc_u::Cdouble = 0.0
  Qhg::Cdouble = 0.0

  # Result part
  Trans::Cdouble = 0.0
  Evap::Cdouble = 0.0
  SH::Cdouble = 0.0
  LH::Cdouble = 0.0
end


function update_ET!(x::OutputET, mid_res::Results, Ta)
  Lv_liquid = (2.501 - 0.00237 * Ta) * 1000000  # The latent heat of water vaporization in j/kg
  
  x.Trans = (x.Trans_o + x.Trans_u) * step

  x.Evap = (x.Eil_o + x.Eil_u +
            x.EiS_o + x.EiS_u +
            x.Evap_soil +
            x.Evap_SW +
            x.Evap_SS) * step
  x.LH = Lv_liquid * (x.Trans_o + x.Trans_u + x.Eil_o + +
                      x.Eil_u + x.Evap_soil + x.Evap_SW) +
         Lv_solid * (x.EiS_o + x.EiS_u + x.Evap_SS)

  x.SH = (x.Qhc_o + x.Qhc_u + x.Qhg)

  # fill values to res
  mid_res.Trans = x.Trans
  mid_res.Evap  = x.Evap
  mid_res.LH    = x.LH
  mid_res.SH    = x.SH
end
