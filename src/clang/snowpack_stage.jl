function snowpack_stage1(Tair, prcp,
  mass_snow_o_last, mass_snow_u_last, mass_snow_g_last,
  mass_snow_o::TypeRef, mass_snow_u::TypeRef, mass_snow_g::TypeRef,
  lai_o::T, lai_u::T, clumping::T,
  area_snow_o::TypeRef, area_snow_u::TypeRef,
  percent_snow_o::TypeRef, percent_snow_u::TypeRef, percent_snow_g::TypeRef,
  density_snow::TypeRef, depth_snow::TypeRef, albedo_v_snow::TypeRef, albedo_n_snow::TypeRef) where {T<:Real}

  # mass_snow_o_last = mass_snow_o[]
  # mass_snow_u_last = mass_snow_u[]
  # mass_snow_g_last = mass_snow_g[]

  ccall((:snowpack_stage1, libbeps), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
    Tair, prcp,
    mass_snow_o_last, mass_snow_u_last, mass_snow_g_last,
    mass_snow_o, mass_snow_u, mass_snow_g, # by reference
    lai_o, lai_u, clumping,
    area_snow_o, area_snow_u, percent_snow_o, percent_snow_u, percent_snow_g, density_snow, depth_snow, albedo_v_snow, albedo_n_snow)

  # percent_snow_o[], percent_snow_u[], percent_snow_g[]
end

function snowpack_stage1(Tair::Float64, prcp::Float64,
  lai_o::Float64, lai_u::Float64, clumping::Float64,
  mass_snow_pre::Layer3{Float64},
  mass_snow::Layer3{Float64}, # by reference
  perc_snow::Layer3{Float64}, # by reference
  area_snow::Layer2{Float64}, # by reference
  depth_snow::Float64, 
  ρ_snow::Ref{Float64}, albedo_v_snow::Ref{Float64}, albedo_n_snow::Ref{Float64})

  _mass_snow_o = Ref(mass_snow.o)
  _mass_snow_u = Ref(mass_snow.u)
  _mass_snow_g = Ref(mass_snow.g)

  _perc_snow_o = Ref(perc_snow.o)
  _perc_snow_u = Ref(perc_snow.u)
  _perc_snow_g = Ref(perc_snow.g)

  _area_snow_o = Ref(area_snow.o)
  _area_snow_u = Ref(area_snow.u)

  _depth_snow = Ref(depth_snow)

  snowpack_stage1(Tair, prcp,
    mass_snow_pre.o, mass_snow_pre.u, mass_snow_pre.g,
    _mass_snow_o, _mass_snow_u, _mass_snow_g,
    lai_o, lai_u, clumping,
    _area_snow_o, _area_snow_u,
    _perc_snow_o, _perc_snow_u, _perc_snow_g,
    ρ_snow, _depth_snow, albedo_v_snow, albedo_n_snow)

  # update values in Layer3
  mass_snow.o = _mass_snow_o[]
  mass_snow.u = _mass_snow_u[]
  mass_snow.g = _mass_snow_g[]
  perc_snow.o = _perc_snow_o[]
  perc_snow.u = _perc_snow_u[]
  perc_snow.g = _perc_snow_g[]
  area_snow.o = _area_snow_o[]
  area_snow.u = _area_snow_u[]
  
  return _depth_snow[]
end



function snowpack_stage2(evapo_snow_o::FT, evapo_snow_u::FT, mass_snow_o::TypeRef, mass_snow_u::TypeRef)
  ccall((:snowpack_stage2, libbeps), Cvoid, (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}),
    evapo_snow_o, evapo_snow_u, mass_snow_o, mass_snow_u)
end

function snowpack_stage2(evapo_snow_o::FT, evapo_snow_u::FT, mass_snow::Layer3{Float64})
  _mass_snow_o = Ref(mass_snow.o)
  _mass_snow_u = Ref(mass_snow.u)
  snowpack_stage2(evapo_snow_o, evapo_snow_u, _mass_snow_o, _mass_snow_u)
  mass_snow.o = _mass_snow_o[]
  mass_snow.u = _mass_snow_u[]
end



function snowpack_stage3(Tair::FT, temp_snow::FT, temp_snow_last::FT, density_snow::FT,
  depth_snow::TypeRef, depth_water::TypeRef, mass_snow_g::TypeRef)
  ccall((:snowpack_stage3, libbeps), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble,
      Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
    Tair, temp_snow, temp_snow_last, density_snow, depth_snow, depth_water, mass_snow_g)
end

function snowpack_stage3(Tair::Float64, Tsnow::Float64, Tsnow_last::Float64, ρ_snow::Float64,
  depth_snow::FT, depth_water::FT, mass_snow::Layer3{Float64})
  # 与julia版本的保持一致
  _depth_snow = Ref(depth_snow)
  _depth_water = Ref(depth_water)
  _mass_snow_g = Ref(mass_snow.g)

  snowpack_stage3(Tair, Tsnow, Tsnow_last, ρ_snow, _depth_snow, _depth_water, _mass_snow_g)
  mass_snow.g = _mass_snow_g[]

  _depth_snow[], _depth_water[]
end


export snowpack_stage3
