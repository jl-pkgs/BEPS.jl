module clang
# using BEPS.beps_c
# using BEPS_jll
# export BEPS_jll
import Base: RefValue, RefArray
TypeRef{T} = Union{Ref{T},RefArray{T,Vector{T},Nothing}}

using BEPS
import BEPS: libbeps

using Parameters
# const libbeps = "deps/libbeps.dll"

# include("../DataType/Constant.jl")
include("struct_SOIL.jl")
include("SOIL_c.jl")
include("module.jl")
# include("beps_inter_prg.jl")
# include("beps_main.jl")
# include("debug_Rln.jl")

# function plantresp(LC, mid_res::Results, lai_yr, lai, Tair, temp_soil, CosZs)
#     ccall((:plantresp, libbeps), Cvoid, (Cint, Ptr{Results}, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble),
#         LC, Ref(mid_res), lai_yr, lai, Tair, temp_soil, CosZs)
# end

function inter_prg_c(jday, rstep,
  lai::T, Ω::T, parameter::Vector{T}, meteo::Met, CosZs::T,
  var_o::Vector{T}, var_n::Vector{T}, soilp::Soil_c,
  Ra::Radiation,
  mid_res::Results, mid_ET::OutputET, var::TransientCache; debug=false, kw...) where {T<:Real}

  ccall((:inter_prg_c, libbeps), Cvoid,
    (Cint, Cint, Cdouble, Cdouble, Ptr{Cdouble},
      Ptr{Met}, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Soil_c}, Ptr{Results}, Ptr{OutputET}),
    jday, rstep, lai, Ω, parameter,
    Ref(meteo), CosZs, var_o, var_n, Ref(soilp), Ref(mid_res), Ref(mid_ET))
end


function Vcmax_Jmax(lai_o, Ω, Vcmax0, slope_Vcmax_N, leaf_N, CosZs, Vcmax_sunlit, Vcmax_shaded, Jmax_sunlit, Jmax_shaded)
  ccall((:Vcmax_Jmax, libbeps), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble,
      Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
    lai_o, Ω, Vcmax0, slope_Vcmax_N, leaf_N, CosZs,
    Vcmax_sunlit, Vcmax_shaded, Jmax_sunlit, Jmax_shaded)
end


# function soilresp(Ccd, Cssd, Csmd, Cfsd, Cfmd, Csm, Cm, Cs, Cp, npp_yr, coef, SoilType, soilp, mid_res)
#     ccall((:soilresp, libbeps), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cfloat, Ptr{Cdouble}, Cint, Ptr{Soil_c}, Ptr{Results}), Ccd, Cssd, Csmd, Cfsd, Cfmd, Csm, Cm, Cs, Cp, npp_yr, coef, SoilType, soilp, mid_res)
# end
function ReadParamVeg(lc::Int=1)
  parameter1 = zeros(Cdouble, 48)
  # ReadParamVeg(short lc, double* parameter1)
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
function lai2(Ω, CosZs, stem_o, stem_u, lai_o, lai_u)
  LAI = Leaf()
  PAI = Leaf()
  lai2(Ω, CosZs, stem_o, stem_u, lai_o, lai_u, LAI, PAI)
  LAI, PAI
end

## add Leaf struct
function lai2(Ω, CosZs, stem_o, stem_u, lai_o, lai_u, LAI::Leaf, PAI::Leaf)
  ccall((:lai2, libbeps), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Leaf}, Ptr{Leaf}),
    Ω, CosZs, stem_o, stem_u, lai_o, lai_u, Ref(LAI), Ref(PAI))
end

function meteo_pack(temp, rh, meteo_pack_output)
  ccall((:meteo_pack, libbeps), Cvoid, (Cdouble, Cdouble, Ptr{Cdouble}), temp, rh, meteo_pack_output)
end

function rainfall_stage1(Tair, prcp, m_water_o_last, m_water_u_last, lai_o, lai_u, Ω)
  m_water_o = init_dbl()
  m_water_u = init_dbl()
  perc_water_o = init_dbl()
  perc_water_u = init_dbl()
  prcp_g = init_dbl()

  ccall((:rainfall_stage1, libbeps), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), Tair, prcp, m_water_o_last, m_water_u_last, lai_o, lai_u, Ω, m_water_o, m_water_u, perc_water_o, perc_water_u, prcp_g)

  m_water_o[], m_water_u[], perc_water_o[], perc_water_u[], prcp_g[]
end

function rainfall_stage1(Tair, prcp, perc_water, m_water, m_water_pre, lai_o, lai_u, Ω)
  # m_water_o_last, m_water_u_last, 
  # perc_water = Layer2{Float64}()
  # m_water = Layer2{Float64}()
  prcp_g = rainfall_stage1_jl(Tair, prcp, perc_water, m_water, m_water_pre, lai_o, lai_u, Ω)
end

function rainfall_stage2(evapo_water_o, evapo_water_u, m_water_o::Ref, m_water_u::Ref)
  ccall((:rainfall_stage2, libbeps), Cvoid, (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), evapo_water_o, evapo_water_u, m_water_o, m_water_u)
end

include("snowpack_stage.jl")

# # no prototype is found for this function at beps.h:181:6, please use with caution
# function rainfall_stage3()
#   ccall((:rainfall_stage3, libbeps), Cvoid, ())
# end

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
  # ReadParamVeg, 
  # readcoef, 
  # lai2, s_coszs,
  # latent_heat!, 
  # Leaf_Temperatures,
  # evaporation_canopy, 
  rainfall_stage1, rainfall_stage2, snowpack_stage1, snowpack_stage2, snowpack_stage3
export Soil_c

end # module
