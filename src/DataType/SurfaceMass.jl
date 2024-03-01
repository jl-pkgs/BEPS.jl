export SurfaceMass


@with_kw mutable struct SurfaceMass{FT<:AbstractFloat}
  ρ_snow       ::FT = 0.0
  prcp_g       ::FT = 0.0

  depth_water  ::FT = 0.0
  depth_snow   ::FT = 0.0

  perc_snow_o  ::FT = 0.0 # Xcs_o
  perc_snow_u  ::FT = 0.0 # Xcs_u
  perc_snow_g::FT = 0.0   # Wcs_g

  perc_water_o ::FT = 0.0 # Xcl_o
  perc_water_u ::FT = 0.0 # Xcl_u

  area_snow_o  ::FT = 0.0 # Ac_snow_o
  area_snow_u  ::FT = 0.0 # Ac_snow_u

  pre_mw_o::FT = 0.0 # Wcl_o
  pre_mw_u::FT = 0.0 # Wcl_u
  
  mw_o::FT = 0.0 # Wcl_o
  mw_u::FT = 0.0 # Wcl_u
  mw_g::FT = 0.0 # Wcl_g

  ms_o  ::FT = 0.0 # Wcs_o
  ms_u  ::FT = 0.0 # Wcs_u
  ms_g  ::FT = 0.0 # Wcs_g

  albedo_v_snow::FT = 0.0
  albedo_n_snow::FT = 0.0
end
