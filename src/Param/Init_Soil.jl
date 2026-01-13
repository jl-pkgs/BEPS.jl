const SOIL_PARAMS = [
  (name="sand",
    b=[1.7, 1.9, 2.1, 2.3, 2.5],
    K_sat=[0.000058, 0.000052, 0.000046, 0.000035, 0.000010],
    θ_sat=0.437, θ_vfc=0.09, θ_vwp=0.03,
    ψ_sat=[0.07, 0.08, 0.09, 0.10, 0.12],
    κ_dry=8.6),
  (name="loamy_sand",
    b=[2.1, 2.3, 2.5, 2.7, 2.9],
    K_sat=[0.000017, 0.000015, 0.000014, 0.000010, 0.000003],
    θ_sat=0.437, θ_vfc=0.21, θ_vwp=0.06,
    ψ_sat=[0.09, 0.10, 0.11, 0.12, 0.14],
    κ_dry=8.3),
  (name="sandy_loam",
    b=[3.1, 3.3, 3.5, 3.7, 3.9],
    K_sat=[0.0000072, 0.00000648, 0.00000576, 0.00000432, 0.00000144],
    θ_sat=0.453, θ_vfc=0.21, θ_vwp=0.10,
    ψ_sat=[0.15, 0.16, 0.17, 0.18, 0.20],
    κ_dry=8.0),
  (name="loam",
    b=[4.5, 4.7, 4.9, 5.1, 5.3],
    K_sat=[0.0000037, 0.0000033, 0.00000296, 0.00000222, 0.00000074],
    θ_sat=0.463, θ_vfc=0.27, θ_vwp=0.12,
    ψ_sat=[0.11, 0.12, 0.13, 0.14, 0.16],
    κ_dry=7.0),
  (name="silty_loam",
    b=[4.7, 4.9, 5.1, 5.3, 5.5],
    K_sat=[0.0000019, 0.0000017, 0.00000152, 0.00000114, 0.00000038],
    θ_sat=0.501, θ_vfc=0.33, θ_vwp=0.13,
    ψ_sat=[0.21, 0.22, 0.23, 0.24, 0.26],
    κ_dry=6.3),
  (name="sandy_clay_loam",
    b=[4.0, 4.2, 4.4, 4.6, 4.8],
    K_sat=[0.0000012, 0.00000108, 0.00000096, 0.00000072, 0.00000024],
    θ_sat=0.398, θ_vfc=0.26, θ_vwp=0.15,
    ψ_sat=[0.28, 0.29, 0.30, 0.31, 0.33],
    κ_dry=7.0),
  (name="clay_loam", # 7
    b=[5.2, 5.4, 5.6, 5.8, 6.0],
    K_sat=[0.00000064, 0.00000058, 0.00000051, 0.00000038, 0.00000013],
    θ_sat=0.464, θ_vfc=0.32, θ_vwp=0.20,
    ψ_sat=[0.26, 0.27, 0.28, 0.29, 0.31],
    κ_dry=5.8),
  (name="silty_clay_loam",
    b=[6.6, 6.8, 7.0, 7.2, 7.4],
    K_sat=[0.00000042, 0.00000038, 0.00000034, 0.000000252, 0.000000084],
    θ_sat=0.471, θ_vfc=0.37, θ_vwp=0.32,
    ψ_sat=[0.33, 0.34, 0.35, 0.36, 0.38],
    κ_dry=4.2),
  (name="sandy_clay",
    b=[6.0, 6.2, 6.4, 6.6, 6.8],
    K_sat=[0.00000033, 0.0000003, 0.000000264, 0.000000198, 0.000000066],
    θ_sat=0.430, θ_vfc=0.34, θ_vwp=0.24,
    ψ_sat=[0.29, 0.30, 0.31, 0.32, 0.34],
    κ_dry=6.3),
  (name="silty_clay",
    b=[7.9, 8.1, 8.3, 8.5, 8.7],
    K_sat=[0.00000025, 0.000000225, 0.0000002, 0.00000015, 0.00000005],
    θ_sat=0.479, θ_vfc=0.39, θ_vwp=0.25,
    ψ_sat=[0.34, 0.35, 0.36, 0.37, 0.39],
    κ_dry=4.0),
  (name="clay",
    b=[7.6, 7.8, 8.0, 8.2, 8.4],
    K_sat=[0.00000017, 0.000000153, 0.000000136, 0.000000102, 0.000000034],
    θ_sat=0.475, θ_vfc=0.40, θ_vwp=0.27,
    ψ_sat=[0.37, 0.38, 0.39, 0.40, 0.42],
    κ_dry=4.4)
]



function init_soil!(soil::Soil, model::BEPSmodel{FT}) where {FT}
  (; hydraulic, thermal, N) = model
  soil.n_layer = Cint(N)

  soil.r_drainage = Cdouble(model.r_drainage)
  soil.r_root_decay = Cdouble(model.r_root_decay)
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


export InitParam_Soil, init_soil!
export SOIL_PARAMS
