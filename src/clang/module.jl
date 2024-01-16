export photosynthesis,
  netRadiation,
  netRadiation_jl,
  evaporation_soil,
  aerodynamic_conductance_c,
  surface_temperature_c,
  aerodynamic_conductance_c,
  transpiration_c,
  sensible_heat_c

# include("s_coszs.jl")
include("evaporation_soil.jl")
include("photosynthesis.jl")
# include("surface_temperature.jl")
# include("netRadiation.jl")
# include("sensible_heat.jl")
# include("transpiration.jl")

function netRadiation(shortRad_global, CosZs,
  temp_o, temp_u, temp_g,
  lai_o, lai_u, lai_os, lai_us, lai::Leaf, Ω, temp_air, rh,
  α_snow_v, α_snow_n,
  percentArea_snow_o, percentArea_snow_u, percent_snow_g,
  α_v_o, α_n_o, α_v_u, α_n_u, α_v_g, α_n_g,
  # netRad_o::TypeRef, netRad_u::TypeRef, netRad_g::TypeRef,
  Rn_Leaf::Leaf, Rns_Leaf::Leaf, Rnl_Leaf::Leaf, Ra::Radiation)

  netRad_o = init_dbl()
  netRad_u = init_dbl()
  netRad_g = init_dbl()

  ccall((:netRadiation, libbeps), Cvoid,
    (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Leaf, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble,
      Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Leaf}, Ptr{Leaf}),
    shortRad_global, CosZs, temp_o, temp_u, temp_g, lai_o, lai_u, lai_os, lai_us, lai, Ω, temp_air, rh, α_snow_v, α_snow_n, percentArea_snow_o, percentArea_snow_u, percent_snow_g,
    α_v_o, α_n_o, α_v_u, α_n_u, α_v_g, α_n_g,
    netRad_o, netRad_u, netRad_g, Ref(Rn_Leaf), Ref(Rns_Leaf))

  netRad_o[], netRad_u[], netRad_g[]
end


function readparam(lc::Int=1)
  parameter1 = zeros(Cdouble, 48)
  # readparam(short lc, double* parameter1)
  ccall((:readparam, libbeps), Cvoid, (Cshort, Ptr{Cdouble}), lc, parameter1)
  parameter1
end

function readcoef(lc::Int=1, stxt::Int=1)
  coef = zeros(Cdouble, 48)
  ccall((:readcoef, libbeps), Cvoid, (Cshort, Cint, Ptr{Cdouble}), lc, stxt, coef)
  coef
end

function aerodynamic_conductance_c(canopy_height_o::T, canopy_height_u::T,
  z_wind::T, clumping::T,
  temp_air::T, wind_sp::T, SH_o_p::T, lai_o::T, lai_u::T=0.0) where {T<:Real}

  ra_o = Ref(0.0)
  ra_u = Ref(0.0)
  ra_g = Ref(0.0)
  G_o_a = Ref(0.0)
  G_o_b = Ref(0.0)
  G_u_a = Ref(0.0)
  G_u_b = Ref(0.0)

  ccall((:aerodynamic_conductance, libbeps), Cvoid,
    (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
    canopy_height_o, canopy_height_u, z_wind, clumping, temp_air, wind_sp, SH_o_p, lai_o, lai_u,
    ra_o, ra_u, ra_g, G_o_a, G_o_b, G_u_a, G_u_b)
  ra_o[], ra_u[], ra_g[], G_o_a[], G_o_b[], G_u_a[], G_u_b[]
end


function transpiration_c(T_leaf::Leaf, Ta::Float64, RH::Float64, Gtrans::Leaf, lai::Leaf)
  trans_o = init_dbl()
  trans_u = init_dbl()
  ccall((:transpiration, libbeps), Cvoid,
    (Leaf, Cdouble, Cdouble, Leaf, Leaf, Ptr{Cdouble}, Ptr{Cdouble}),
    T_leaf, Ta, RH, Gtrans, lai, trans_o, trans_u)
  trans_o[], trans_u[]
end

function sensible_heat_c(tempL::Leaf,
  temp_g::Cdouble, temp_air::Cdouble, RH::Cdouble,
  Gheat::Leaf, Gheat_g::Cdouble, LAI::Leaf)

  SH_o = init_dbl()
  SH_u = init_dbl()
  SH_g = init_dbl()

  ccall((:sensible_heat, libbeps), Cvoid,
    (Leaf, Cdouble, Cdouble, Cdouble, Leaf, Cdouble, Leaf,
      Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
    tempL, temp_g, temp_air, RH, Gheat, Gheat_g, LAI,
    SH_o, SH_u, SH_g)
  SH_o[], SH_u[], SH_g[]
end

function surface_temperature_c(T_air::FT, rh_air::FT, depth_snow::FT, depth_water::FT,
  capacity_heat_soil1::FT, capacity_heat_soil0::FT, Gheat_g::FT,
  depth_soil1::FT, density_snow::FT, tempL_u::FT, netRad_g::FT,
  evapo_soil::FT, evapo_water_g::FT, evapo_snow_g::FT, lambda_soil1::FT,
  percent_snow_g::FT, heat_flux_soil1::FT, T_ground_last::FT,
  T_soil1_last::FT, T_any0_last::FT, T_soil0_last::FT,
  T_snow_last::FT, T_snow1_last::FT, T_snow2_last::FT)

  T_ground = Ref(0.0)
  T_any0 = Ref(0.0)
  T_snow = Ref(0.0)
  T_soil0 = Ref(0.0)
  T_snow1 = Ref(0.0)
  T_snow2 = Ref(0.0)
  heat_flux = Ref(0.0)

  ccall((:surface_temperature, libbeps), Cvoid,
    (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble,
      Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
    T_air, rh_air, depth_snow, depth_water, capacity_heat_soil1, capacity_heat_soil0, Gheat_g, depth_soil1, density_snow, tempL_u, netRad_g, evapo_soil, evapo_water_g, evapo_snow_g, lambda_soil1, percent_snow_g, heat_flux_soil1, T_ground_last, T_soil1_last, T_any0_last, T_snow_last, T_soil0_last, T_snow1_last, T_snow2_last, T_ground, T_any0, T_snow, T_soil0, T_snow1, T_snow2, heat_flux)

  heat_flux[], T_ground[], T_any0[], T_soil0[], T_snow[], T_snow1[], T_snow2[]
end
