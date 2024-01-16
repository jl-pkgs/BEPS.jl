module clang
# using BEPS.beps_c
# using BEPS_jll
# export BEPS_jll
import Base: RefValue, RefArray
TypeRef{T} = Union{Ref{T},RefArray{T,Vector{T},Nothing}}

using BEPS
import BEPS: libbeps
# const libbeps = "deps/libbeps.dll"

# include("../DataType/Constant.jl")
include("SOIL_c.jl")
include("module.jl")
# include("beps_inter_prg.jl")
# include("beps_main.jl")
# include("debug_Rln.jl")

# function plantresp(LC, mid_res::Results, lai_yr, lai, temp_air, temp_soil, CosZs)
#     ccall((:plantresp, libbeps), Cvoid, (Cint, Ptr{Results}, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble),
#         LC, Ref(mid_res), lai_yr, lai, temp_air, temp_soil, CosZs)
# end

function inter_prg_c(jday, rstep,
  lai::T, clumping::T, parameter::Vector{T}, meteo::ClimateData, CosZs::T,
  var_o::Vector{T}, var_n::Vector{T}, soilp::Soil_c,
  Ra::Radiation,
  mid_res::Results, mid_ET::OutputET, var::InterTempVars; debug=false, kw...) where {T<:Real}

  ccall((:inter_prg_c, libbeps), Cvoid,
    (Cint, Cint, Cdouble, Cdouble, Ptr{Cdouble},
      Ptr{ClimateData}, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Soil_c}, Ptr{Results}, Ptr{OutputET}),
    jday, rstep, lai, clumping, parameter,
    Ref(meteo), CosZs, var_o, var_n, Ref(soilp), Ref(mid_res), Ref(mid_ET))
end


function Vcmax_Jmax(lai_o, clumping, Vcmax0, slope_Vcmax_N, leaf_N, CosZs, Vcmax_sunlit, Vcmax_shaded, Jmax_sunlit, Jmax_shaded)
  ccall((:Vcmax_Jmax, libbeps), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble,
      Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
    lai_o, clumping, Vcmax0, slope_Vcmax_N, leaf_N, CosZs,
    Vcmax_sunlit, Vcmax_shaded, Jmax_sunlit, Jmax_shaded)
end


# function soilresp(Ccd, Cssd, Csmd, Cfsd, Cfmd, Csm, Cm, Cs, Cp, npp_yr, coef, soiltype, soilp, mid_res)
#     ccall((:soilresp, libbeps), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cfloat, Ptr{Cdouble}, Cint, Ptr{Soil_c}, Ptr{Results}), Ccd, Cssd, Csmd, Cfsd, Cfmd, Csm, Cm, Cs, Cp, npp_yr, coef, soiltype, soilp, mid_res)
# end
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

"""
s_coszs(jday::Int, j::Int, lat::T, lon::T) where {T<:Real}

# Example
```jldoc
x = zeros(4)
s_coszs(Ref(x, 1), 10, 1, 30, 120)
x
```
"""
function s_coszs(jday::Int, j::Int, lat::T, lon::T) where {T<:Real}
  CosZs = init_dbl()
  ccall((:s_coszs, libbeps), Cvoid, (Cshort, Cshort, Cfloat, Cfloat, Ptr{Cdouble}),
    jday, j, lat, lon, CosZs)
  # return Value(CosZs)
  CosZs[]
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

function meteo_pack(temp, rh, meteo_pack_output)
  ccall((:meteo_pack, libbeps), Cvoid, (Cdouble, Cdouble, Ptr{Cdouble}), temp, rh, meteo_pack_output)
end

function rainfall_stage1(temp_air, prcp, mass_water_o_last, mass_water_u_last, lai_o, lai_u, clumping)
  mass_water_o = init_dbl()
  mass_water_u = init_dbl()
  percent_water_o = init_dbl()
  percent_water_u = init_dbl()
  prcp_g = init_dbl()

  ccall((:rainfall_stage1, libbeps), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), temp_air, prcp, mass_water_o_last, mass_water_u_last, lai_o, lai_u, clumping, mass_water_o, mass_water_u, percent_water_o, percent_water_u, prcp_g)

  mass_water_o[], mass_water_u[], percent_water_o[], percent_water_u[], prcp_g[]
end

function rainfall_stage2(evapo_water_o, evapo_water_u, mass_water_o::Ref, mass_water_u::Ref)
  ccall((:rainfall_stage2, libbeps), Cvoid, (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), evapo_water_o, evapo_water_u, mass_water_o, mass_water_u)
end

# # no prototype is found for this function at beps.h:181:6, please use with caution
# function rainfall_stage3()
#   ccall((:rainfall_stage3, libbeps), Cvoid, ())
# end

function snowpack_stage1(temp_air, prcp,
  # mass_snow_o_last, mass_snow_u_last, mass_snow_g_last,
  mass_snow_o::TypeRef, mass_snow_u::TypeRef, mass_snow_g::TypeRef,
  lai_o::T, lai_u::T, clumping::T,
  area_snow_o::TypeRef, area_snow_u::TypeRef, 
  # percent_snow_o::TypeRef, percent_snow_u::TypeRef, percent_snow_g::TypeRef,
  density_snow::TypeRef, depth_snow::TypeRef, albedo_v_snow::TypeRef, albedo_n_snow::TypeRef) where {T<:Real}

  mass_snow_o_last = mass_snow_o[]
  mass_snow_u_last = mass_snow_u[]
  mass_snow_g_last = mass_snow_g[]
  percent_snow_o = init_dbl()
  percent_snow_u = init_dbl()
  percent_snow_g = init_dbl()

  ccall((:snowpack_stage1, libbeps), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
    temp_air, prcp,
    mass_snow_o_last, mass_snow_u_last, mass_snow_g_last,
    mass_snow_o, mass_snow_u, mass_snow_g, # by reference
    lai_o, lai_u, clumping,
    area_snow_o, area_snow_u, percent_snow_o, percent_snow_u, percent_snow_g, density_snow, depth_snow, albedo_v_snow, albedo_n_snow)

  percent_snow_o[], percent_snow_u[], percent_snow_g[]
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
# for name in names(@__MODULE__; all=true)
#   # if startswith(string(name), prefix)
#   if !(string(name) in ["#eval", "#include"])
#     @eval export $name
#   end
#   # end
# end

export inter_prg_c,
  # Vcmax_Jmax, 
  # readparam, 
  # readcoef, 
  # lai2, s_coszs,
  # latent_heat!, 
  # Leaf_Temperatures,
  # evaporation_canopy, 
  rainfall_stage1, rainfall_stage2, snowpack_stage1, snowpack_stage2, snowpack_stage3


end # module
