module beps_c
# using BEPS.beps_c
# using BEPS_jll
# export BEPS_jll

# using CEnum
const libbeps = "deps/libbeps.dll"


import Base: RefValue, RefArray

TypeRef{T} = Union{Ref{T},RefArray{T,Vector{T},Nothing}}

include("DataType/DataType.jl")
include("helpers.jl")
include("SOIL_c.jl")
include("module/module.jl")
include("beps_inter_prg.jl")
include("beps_main.jl")

# function plantresp(LC, mid_res::Results, lai_yr, lai, temp_air, temp_soil, CosZs)
#     ccall((:plantresp, libbeps), Cvoid, (Cint, Ptr{Results}, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble),
#         LC, Ref(mid_res), lai_yr, lai, temp_air, temp_soil, CosZs)
# end

function Vcmax_Jmax(lai_o, clumping, Vcmax0, slope_Vcmax_N, leaf_N, CosZs, Vcmax_sunlit, Vcmax_shaded, Jmax_sunlit, Jmax_shaded)
  ccall((:Vcmax_Jmax, libbeps), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble,
      Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
    lai_o, clumping, Vcmax0, slope_Vcmax_N, leaf_N, CosZs,
    Vcmax_sunlit, Vcmax_shaded, Jmax_sunlit, Jmax_shaded)
end


# function soilresp(Ccd, Cssd, Csmd, Cfsd, Cfmd, Csm, Cm, Cs, Cp, npp_yr, coef, soiltype, soilp, mid_res)
#     ccall((:soilresp, libbeps), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cfloat, Ptr{Cdouble}, Cint, Ptr{Soil}, Ptr{Results}), Ccd, Cssd, Csmd, Cfsd, Cfmd, Csm, Cm, Cs, Cp, npp_yr, coef, soiltype, soilp, mid_res)
# end

# √, passed test
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


# no prototype is found for this function at beps.h:133:6, please use with caution
function readhydr_param()
  ccall((:readhydr_param, libbeps), Cvoid, ())
end



# no prototype is found for this function at beps.h:138:6, please use with caution
function soil_water_factor()
  ccall((:soil_water_factor, libbeps), Cvoid, ())
end

# passed test
# lai2(0.8, 0.6, 0.2, 0.4, 0.2, 0.4)
function lai2(clumping, CosZs, stem_o, stem_u, lai_o, lai_u)
  LAI = Leaf()
  PAI = Leaf()
  lai2(clumping, CosZs, stem_o, stem_u, lai_o, lai_u, LAI, PAI)
  LAI, PAI
end

## add Leaf struct
function lai2(clumping, CosZs, stem_o, stem_u, lai_o, lai_u, LAI::Leaf, PAI::Leaf)
  ccall((:lai2, libbeps), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Leaf}, Ptr{Leaf}),
    clumping, CosZs, stem_o, stem_u, lai_o, lai_u, Ref(LAI), Ref(PAI))
end

function Leaf_Temperatures(Tair, slope, psychrometer, VPD_air, Cp_ca, Gw, Gww, Gh, Xcs_o, Xcl_o, Xcs_u, Xcl_u, radiation::Leaf, Tc::Leaf)
  ccall((:Leaf_Temperatures, libbeps), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Leaf, Leaf, Leaf, Cdouble, Cdouble, Cdouble, Cdouble, Leaf, Ptr{Leaf}),
    Tair, slope, psychrometer, VPD_air, Cp_ca, Gw, Gww, Gh, 
    Xcs_o, Xcl_o, Xcs_u, Xcl_u, 
    radiation, Ref(Tc))
end

function Leaf_Temperature(Tair, slope, psychrometer, VPD_air, Cp_ca, Gw, Gww, Gh, Xcs, Xcl, radiation)
  ccall((:Leaf_Temperature, libbeps), Cdouble, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble), 
  Tair, slope, psychrometer, VPD_air, Cp_ca, Gw, Gww, Gh, Xcs, Xcl, radiation)
end


function latent_heat!(leleaf::Leaf, Gw::Leaf, VPD_air, slope, Tc_old::Leaf, temp_air, rho_a, Cp_ca, psychrometer)
  leleaf.o_sunlit = Gw.o_sunlit * (VPD_air + slope * (Tc_old.o_sunlit - temp_air)) * rho_a * Cp_ca / psychrometer
  leleaf.o_shaded = Gw.o_shaded * (VPD_air + slope * (Tc_old.o_shaded - temp_air)) * rho_a * Cp_ca / psychrometer
  leleaf.u_sunlit = Gw.u_sunlit * (VPD_air + slope * (Tc_old.u_sunlit - temp_air)) * rho_a * Cp_ca / psychrometer
  leleaf.u_shaded = Gw.u_shaded * (VPD_air + slope * (Tc_old.u_shaded - temp_air)) * rho_a * Cp_ca / psychrometer
end


function evaporation_canopy(tempL, temp_air, rh_air, Gwater, lai, percent_water_o, percent_water_u, percent_snow_o, percent_snow_u, evapo_water_o, evapo_water_u, evapo_snow_o, evapo_snow_u)
  ccall((:evaporation_canopy, libbeps), Cvoid, (Leaf, Cdouble, Cdouble, Leaf, Leaf, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), tempL, temp_air, rh_air, Gwater, lai, percent_water_o, percent_water_u, percent_snow_o, percent_snow_u, evapo_water_o, evapo_water_u, evapo_snow_o, evapo_snow_u)
end

# end of leaf struct


function rainfall_stage1(temp_air, precipitation, mass_water_o_last, mass_water_u_last, lai_o, lai_u, clumping, mass_water_o, mass_water_u, percent_water_o, percent_water_u, precipitation_g)
  ccall((:rainfall_stage1, libbeps), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), temp_air, precipitation, mass_water_o_last, mass_water_u_last, lai_o, lai_u, clumping, mass_water_o, mass_water_u, percent_water_o, percent_water_u, precipitation_g)
end

function rainfall_stage2(evapo_water_o, evapo_water_u, mass_water_o, mass_water_u)
  ccall((:rainfall_stage2, libbeps), Cvoid, (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), evapo_water_o, evapo_water_u, mass_water_o, mass_water_u)
end

# no prototype is found for this function at beps.h:181:6, please use with caution
function rainfall_stage3()
  ccall((:rainfall_stage3, libbeps), Cvoid, ())
end

function meteo_pack(temp, rh, meteo_pack_output)
  ccall((:meteo_pack, libbeps), Cvoid, (Cdouble, Cdouble, Ptr{Cdouble}), temp, rh, meteo_pack_output)
end


function snowpack_stage1(temp_air, precipitation,
  mass_snow_o_last, mass_snow_u_last, mass_snow_g_last,
  mass_snow_o::TypeRef, mass_snow_u::TypeRef, mass_snow_g::TypeRef,
  lai_o::T, lai_u::T, clumping::T,
  area_snow_o::TypeRef, area_snow_u::TypeRef, percent_snow_o::TypeRef, percent_snow_u::TypeRef, percent_snow_g::TypeRef,
  density_snow::TypeRef, depth_snow::TypeRef, albedo_v_snow::TypeRef, albedo_n_snow::TypeRef) where {T<:Real}

  ccall((:snowpack_stage1, libbeps), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
    temp_air, precipitation,
    mass_snow_o_last, mass_snow_u_last, mass_snow_g_last,
    mass_snow_o, mass_snow_u, mass_snow_g, # by reference
    lai_o, lai_u, clumping,
    area_snow_o, area_snow_u, percent_snow_o, percent_snow_u, percent_snow_g, density_snow, depth_snow, albedo_v_snow, albedo_n_snow)
end

function snowpack_stage2(evapo_snow_o, evapo_snow_u, mass_snow_o, mass_snow_u)
  ccall((:snowpack_stage2, libbeps), Cvoid, (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}),
    evapo_snow_o, evapo_snow_u, mass_snow_o, mass_snow_u)
end

function snowpack_stage3(temp_air, temp_snow, temp_snow_last, density_snow, depth_snow, depth_water, mass_snow_g)
  ccall((:snowpack_stage3, libbeps), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
    temp_air, temp_snow, temp_snow_last, density_snow, depth_snow, depth_water, mass_snow_g)
end

# exports
# const PREFIXES = ["CX", "clang_"]
# , prefix in PREFIXES
for name in names(@__MODULE__; all=true)
  # if startswith(string(name), prefix)
  if !(string(name) in ["eval", "include"])
    @eval export $name
  end
  # end
end

end # module
