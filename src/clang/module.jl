export 
  evaporation_soil_c,
  evaporation_canopy_c,
  netRadiation_c,
  photosynthesis_c,
  aerodynamic_conductance_c,
  Leaf_Temperature_c,
  surface_temperature_c,
  aerodynamic_conductance_c,
  transpiration_c,
  sensible_heat_c

# include("s_coszs.jl")
# include("evaporation_soil.jl")
# include("surface_temperature.jl")
# include("netRadiation.jl")
# include("sensible_heat.jl")
# include("transpiration.jl")
function photosynthesis_c(temp_leaf_p::Cdouble, rad_leaf::Cdouble, e_air::Cdouble,
  g_lb_w::Cdouble, vc_opt::Cdouble,
  f_soilwater::Cdouble, b_h2o::Cdouble, m_h2o::Cdouble,
  cii::Cdouble, temp_leaf_c::Cdouble, LH_leaf::Cdouble)

  Gs_w = init_dbl()
  aphoto = init_dbl()
  ci = init_dbl()

  ccall((:photosynthesis, libbeps), Cvoid,
    (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble,
      Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
    temp_leaf_p, rad_leaf, e_air, g_lb_w, vc_opt, f_soilwater, b_h2o, m_h2o, cii, temp_leaf_c, LH_leaf,
    Gs_w, aphoto, ci)

  Gs_w[], aphoto[], ci[]
end

function netRadiation_c(shortRad_global, CosZs,
  temp_o, temp_u, temp_g,
  lai_o, lai_u, lai_os, lai_us, lai::Leaf, Ω, temp_air, rh,
  α_snow_v, α_snow_n,
  percArea_snow_o, percArea_snow_u, perc_snow_g,
  α_v_o, α_n_o, α_v_u, α_n_u, α_v_g, α_n_g,
  # netRad_o::TypeRef, netRad_u::TypeRef, netRad_g::TypeRef,
  Rn_Leaf::Leaf, Rns_Leaf::Leaf, Rnl_Leaf::Leaf, Ra::Radiation)

  netRad_o = init_dbl()
  netRad_u = init_dbl()
  netRad_g = init_dbl()

  ccall((:netRadiation, libbeps), Cvoid,
    (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Leaf, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble,
      Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Leaf}, Ptr{Leaf}),
    shortRad_global, CosZs, temp_o, temp_u, temp_g, lai_o, lai_u, lai_os, lai_us, lai, Ω, temp_air, rh, α_snow_v, α_snow_n, percArea_snow_o, percArea_snow_u, perc_snow_g,
    α_v_o, α_n_o, α_v_u, α_n_u, α_v_g, α_n_g,
    netRad_o, netRad_u, netRad_g, Ref(Rn_Leaf), Ref(Rns_Leaf))

  netRad_o[], netRad_u[], netRad_g[]
end


function evaporation_soil_c(temp_air, temp_g, rh_air, netRad_g, Gheat_g,
  perc_snow_g::TypeRef, z_water::TypeRef, z_snow::TypeRef, mass_water_g::TypeRef, mass_snow_g::TypeRef,
  density_snow, swc_g, porosity_g)
  # evapo_soil::TypeRef, evapo_water_g::TypeRef, evapo_snow_g::TypeRef

  evapo_soil = init_dbl()
  evapo_water_g = init_dbl()
  evapo_snow_g = init_dbl()

  ccall((:evaporation_soil, libbeps), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble, Cdouble,
      Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
    temp_air, temp_g, rh_air, netRad_g, Gheat_g, perc_snow_g, z_water, z_snow, mass_water_g, mass_snow_g, density_snow, swc_g, porosity_g, evapo_soil, evapo_water_g, evapo_snow_g)

  evapo_soil[], evapo_water_g[], evapo_snow_g[]
end

function evaporation_canopy_c(tempL::Leaf, Ta::Float64, rh_air::Float64,
  Gwater::Leaf, lai::Leaf,
  perc_water_o::Float64, perc_water_u::Float64, perc_snow_o::Float64, perc_snow_u::Float64)

  evapo_water_o = init_dbl()
  evapo_water_u = init_dbl()
  evapo_snow_o = init_dbl()
  evapo_snow_u = init_dbl()

  ccall((:evaporation_canopy, libbeps), Cvoid, (Leaf, Cdouble, Cdouble, Leaf, Leaf, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), tempL, Ta, rh_air, Gwater, lai, perc_water_o, perc_water_u, perc_snow_o, perc_snow_u, evapo_water_o, evapo_water_u, evapo_snow_o, evapo_snow_u)

  evapo_water_o[], evapo_water_u[], evapo_snow_o[], evapo_snow_u[]
end

function Leaf_Temperatures_c(Tair, slope, γ, VPD_air, Cp_ca, Gw, Gww, Gh, Xcs_o, Xcl_o, Xcs_u, Xcl_u, radiation::Leaf, Tc::Leaf)
  ccall((:Leaf_Temperatures, libbeps), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Leaf, Leaf, Leaf, Cdouble, Cdouble, Cdouble, Cdouble, Leaf, Ptr{Leaf}),
    Tair, slope, γ, VPD_air, Cp_ca, Gw, Gww, Gh,
    Xcs_o, Xcl_o, Xcs_u, Xcl_u,
    radiation, Ref(Tc))
end

function Leaf_Temperature_c(Tair, slope, γ, VPD_air, Cp_ca, Gw, Gww, Gh, Xcs, Xcl, radiation)
  ccall((:Leaf_Temperature, libbeps), Cdouble, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble),
    Tair, slope, γ, VPD_air, Cp_ca, Gw, Gww, Gh, Xcs, Xcl, radiation)
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

function surface_temperature_c(T_air::FT, rh_air::FT, z_snow::FT, z_water::FT,
  capacity_heat_soil1::FT, capacity_heat_soil0::FT, Gheat_g::FT,
  depth_soil1::FT, density_snow::FT, tempL_u::FT, netRad_g::FT,
  evapo_soil::FT, evapo_water_g::FT, evapo_snow_g::FT, lambda_soil1::FT,
  perc_snow_g::FT, heat_flux_soil1::FT, T_ground_last::FT,
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
    T_air, rh_air, z_snow, z_water, capacity_heat_soil1, capacity_heat_soil0, Gheat_g, depth_soil1, density_snow, tempL_u, netRad_g, evapo_soil, evapo_water_g, evapo_snow_g, lambda_soil1, perc_snow_g, heat_flux_soil1, T_ground_last, T_soil1_last, T_any0_last, T_snow_last, T_soil0_last, T_snow1_last, T_snow2_last, T_ground, T_any0, T_snow, T_soil0, T_snow1, T_snow2, heat_flux)

  heat_flux[], T_ground[], T_any0[], T_soil0[], T_snow[], T_snow1[], T_snow2[]
end
