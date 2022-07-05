module beps_c
# using BEPS.beps_c
# using BEPS_jll
# export BEPS_jll

# using CEnum
const libbeps = "lib/libbeps.dll"

include("SOIL_c.jl")
include("beps_main.jl")

Value = getindex

# init *Cdouble
init_dbl() = Ref(0.0)

function inter_prg(jday, rstep, lai, clumping, parameter,
    meteo::ClimateData, CosZs, var_o, var_n, soilp::Soil, mid_res::Results)

    ccall((:inter_prg, libbeps), Cvoid,
        (Cint, Cint, Cdouble, Cdouble, Ptr{Cdouble},
            Ptr{ClimateData}, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Soil}, Ptr{Results}),
        jday, rstep, lai, clumping, parameter,
        Ref(meteo), CosZs, var_o, var_n, Ref(soilp), Ref(mid_res))
end

function s_coszs(jday::Int, j::Int, lat::T, lon::T) where {T<:Real}
    ## This version works
    # CosZs = [0.0]
    # ccall((:s_coszs, libbeps), Cvoid, (Cshort, Cshort, Cfloat, Cfloat, Ptr{Cdouble}),
    #     jday, j, lat, lon, CosZs)
    CosZs = init_dbl()
    ccall((:s_coszs, libbeps), Cvoid, (Cshort, Cshort, Cfloat, Cfloat, Ptr{Cdouble}),
        jday, j, lat, lon, CosZs)
    # Cshort(jday), Cshort(j), Cfloat(lat), Cfloat(lon), CosZs)
    return Value(CosZs)
end

function aerodynamic_conductance(canopy_height_o, canopy_height_u, zz, clumping, temp_air, wind_sp, SH_o_p, lai_o, lai_u,
    rm, ra_u, ra_g, G_o_a, G_o_b, G_u_a, G_u_b)

    ccall((:aerodynamic_conductance, libbeps), Cvoid,
        (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
        canopy_height_o, canopy_height_u, zz, clumping, temp_air, wind_sp, SH_o_p, lai_o, lai_u,
        rm, ra_u, ra_g, G_o_a, G_o_b, G_u_a, G_u_b)
end

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

function photosynthesis(temp_leaf_p, rad_leaf, e_air, g_lb_w, vc_opt, f_soilwater, b_h2o, m_h2o, cii, temp_leaf_c, LH_leaf, Gs_w, aphoto, ci)
    ccall((:photosynthesis, libbeps), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), temp_leaf_p, rad_leaf, e_air, g_lb_w, vc_opt, f_soilwater, b_h2o, m_h2o, cii, temp_leaf_c, LH_leaf, Gs_w, aphoto, ci)
end

# no prototype is found for this function at beps.h:138:6, please use with caution
function soil_water_factor()
    ccall((:soil_water_factor, libbeps), Cvoid, ())
end

## add Leaf struct
function lai2(clumping, CosZs, stem_o, stem_u, lai_o, lai_u, LAI, PAI)
    ccall((:lai2, libbeps), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Leaf}, Ptr{Leaf}),
        clumping, CosZs, stem_o, stem_u, lai_o, lai_u, LAI, PAI)
end

function netRadiation(shortRad_global, CosZs, temp_o, temp_u, temp_g, lai_o, lai_u, lai_os, lai_us, lai, clumping, temp_air, rh,
    albedo_snow_v, albedo_snow_n,
    percentArea_snow_o, percentArea_snow_u, percent_snow_g,
    albedo_v_o, albedo_n_o, albedo_v_u, albedo_n_u, albedo_v_g, albedo_n_g,
    netRad_o, netRad_u, netRad_g, netRadLeaf, netShortRadLeaf)

    ccall((:netRadiation, libbeps), Cvoid,
        (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Leaf, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Leaf}, Ptr{Leaf}),
        shortRad_global, CosZs, temp_o, temp_u, temp_g, lai_o, lai_u, lai_os, lai_us, lai, clumping, temp_air, rh, albedo_snow_v, albedo_snow_n, percentArea_snow_o, percentArea_snow_u, percent_snow_g, 
        albedo_v_o, albedo_n_o, albedo_v_u, albedo_n_u, albedo_v_g, albedo_n_g,
        netRad_o, netRad_u, netRad_g, netRadLeaf, netShortRadLeaf)
end

function Leaf_Temperatures(Tair, slope, psychrometer, VPD_air, Cp_ca, Gw, Gww, Gh, Xcs_o, Xcl_o, Xcs_u, Xcl_u, radiation, Tc)
    ccall((:Leaf_Temperatures, libbeps), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Leaf, Leaf, Leaf, Cdouble, Cdouble, Cdouble, Cdouble, Leaf, Ptr{Leaf}),
        Tair, slope, psychrometer, VPD_air, Cp_ca, Gw, Gww, Gh, Xcs_o, Xcl_o, Xcs_u, Xcl_u, radiation, Tc)
end

function Leaf_Temperature(Tair, slope, psychrometer, VPD_air, Cp_ca, Gw, Gww, Gh, Xcs, Xcl, radiation)
    ccall((:Leaf_Temperature, libbeps), Cdouble, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble), Tair, slope, psychrometer, VPD_air, Cp_ca, Gw, Gww, Gh, Xcs, Xcl, radiation)
end

function sensible_heat(tempL, temp_g, temp_air, rh_air, Gheat, Gheat_g, LAI, SH_o, SH_u, SH_g)
    ccall((:sensible_heat, libbeps), Cvoid, (Leaf, Cdouble, Cdouble, Cdouble, Leaf, Cdouble, Leaf, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), tempL, temp_g, temp_air, rh_air, Gheat, Gheat_g, LAI, SH_o, SH_u, SH_g)
end

function transpiration(tempL, temp_air, rh_air, Gtrans, lai, trans_o, trans_u)
    ccall((:transpiration, libbeps), Cvoid, (Leaf, Cdouble, Cdouble, Leaf, Leaf, Ptr{Cdouble}, Ptr{Cdouble}), tempL, temp_air, rh_air, Gtrans, lai, trans_o, trans_u)
end

function evaporation_canopy(tempL, temp_air, rh_air, Gwater, lai, percent_water_o, percent_water_u, percent_snow_o, percent_snow_u, evapo_water_o, evapo_water_u, evapo_snow_o, evapo_snow_u)
    ccall((:evaporation_canopy, libbeps), Cvoid, (Leaf, Cdouble, Cdouble, Leaf, Leaf, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), tempL, temp_air, rh_air, Gwater, lai, percent_water_o, percent_water_u, percent_snow_o, percent_snow_u, evapo_water_o, evapo_water_u, evapo_snow_o, evapo_snow_u)
end

# end of leaf struct

function evaporation_soil(temp_air, temp_g, rh_air, netRad_g, Gheat_g, percent_snow_g, depth_water, depth_snow, mass_water_g, mass_snow_g, density_snow, swc_g, porosity_g, evapo_soil, evapo_water_g, evapo_snow_g)
    ccall((:evaporation_soil, libbeps), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), temp_air, temp_g, rh_air, netRad_g, Gheat_g, percent_snow_g, depth_water, depth_snow, mass_water_g, mass_snow_g, density_snow, swc_g, porosity_g, evapo_soil, evapo_water_g, evapo_snow_g)
end

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

function surface_temperature(temp_air, rh_air, depth_snow, depth_water, capacity_heat_soil1, capacity_heat_soil0, Gheat_g, depth_soil1, density_snow, tempL_u, netRad_g, evapo_soil, evapo_water_g, evapo_snow_g, lambda_soil1, percent_snow_g, heat_flux_soil1, temp_ground_last, temp_soil1_last, temp_any0_last, temp_snow_last, temp_soil0_last, temp_snow1_last, temp_snow2_last, temp_ground, temp_any0, temp_snow, temp_soil0, temp_snow1, temp_snow2, heat_flux)
    ccall((:surface_temperature, libbeps), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), temp_air, rh_air, depth_snow, depth_water, capacity_heat_soil1, capacity_heat_soil0, Gheat_g, depth_soil1, density_snow, tempL_u, netRad_g, evapo_soil, evapo_water_g, evapo_snow_g, lambda_soil1, percent_snow_g, heat_flux_soil1, temp_ground_last, temp_soil1_last, temp_any0_last, temp_snow_last, temp_soil0_last, temp_snow1_last, temp_snow2_last, temp_ground, temp_any0, temp_snow, temp_soil0, temp_snow1, temp_snow2, heat_flux)
end

function snowpack_stage1(temp_air, precipitation, mass_snow_o_last, mass_snow_u_last, mass_snow_g_last, mass_snow_o, mass_snow_u, mass_snow_g, lai_o, lai_u, clumping, area_snow_o, area_snow_u, percent_snow_o, percent_snow_u, percent_snow_g, density_snow, depth_snow, albedo_v_snow, albedo_n_snow)
    ccall((:snowpack_stage1, libbeps), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
        temp_air, precipitation, mass_snow_o_last, mass_snow_u_last, mass_snow_g_last, mass_snow_o, mass_snow_u, mass_snow_g, lai_o, lai_u, clumping, area_snow_o, area_snow_u, percent_snow_o, percent_snow_u, percent_snow_g, density_snow, depth_snow, albedo_v_snow, albedo_n_snow)
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
