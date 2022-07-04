module BEPS

using BEPS_jll
export BEPS_jll

using CEnum

include("SOIL_c.jl")

function inter_prg(jday, rstep, lai, clumping, parameter, meteo, CosZs, var_o, var_n, soilp, mid_res)
    ccall((:inter_prg, libbeps), Cvoid, (Cint, Cint, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{ClimateData}, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Soil}, Ptr{Results}), 
    jday, rstep, lai, clumping, parameter, meteo, CosZs, var_o, var_n, soilp, mid_res)
end

function s_coszs(jday, j, lat, lon, CosZs)
    ccall((:s_coszs, libbeps), Cvoid, (Cshort, Cshort, Cfloat, Cfloat, Ptr{Cdouble}), 
    jday, j, lat, lon, CosZs)
end

function aerodynamic_conductance(canopy_height_o, canopy_height_u, zz, clumping, temp_air, wind_sp, SH_o_p, lai_o, lai_u, rm, ra_u, ra_g, G_o_a, G_o_b, G_u_a, G_u_b)
    ccall((:aerodynamic_conductance, libbeps), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), canopy_height_o, canopy_height_u, zz, clumping, temp_air, wind_sp, SH_o_p, lai_o, lai_u, rm, ra_u, ra_g, G_o_a, G_o_b, G_u_a, G_u_b)
end

function plantresp(LC, mid_res, lai_yr, lai, temp_air, temp_soil, CosZs)
    ccall((:plantresp, libbeps), Cvoid, (Cint, Ptr{Results}, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble), LC, mid_res, lai_yr, lai, temp_air, temp_soil, CosZs)
end

function Vcmax_Jmax(lai_o, clumping, Vcmax0, slope_Vcmax_N, leaf_N, CosZs, Vcmax_sunlit, Vcmax_shaded, Jmax_sunlit, Jmax_shaded)
    ccall((:Vcmax_Jmax, libbeps), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), lai_o, clumping, Vcmax0, slope_Vcmax_N, leaf_N, CosZs, Vcmax_sunlit, Vcmax_shaded, Jmax_sunlit, Jmax_shaded)
end

function netRadiation(shortRad_global, CosZs, temp_o, temp_u, temp_g, lai_o, lai_u, lai_os, lai_us, lai_o_sunlit, lai_o_shaded, lai_u_sunlit, lai_u_shaded, clumping, temp_air, rh, albedo_snow_v, albedo_snow_n, percentArea_snow_o, percentArea_snow_u, percent_snow_g, albedo_v_o, albedo_n_o, albedo_v_u, albedo_n_u, albedo_v_g, albedo_n_g, netRad_o, netRad_u, netRad_g, netRadLeaf_o_sunlit, netRadLeaf_o_shaded, netRadLeaf_u_sunlit, netRadLeaf_u_shaded, netShortRadLeaf_o_sunlit, netShortRadLeaf_o_shaded, netShortRadLeaf_u_sunlit, netShortRadLeaf_u_shaded)
    ccall((:netRadiation, libbeps), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), shortRad_global, CosZs, temp_o, temp_u, temp_g, lai_o, lai_u, lai_os, lai_us, lai_o_sunlit, lai_o_shaded, lai_u_sunlit, lai_u_shaded, clumping, temp_air, rh, albedo_snow_v, albedo_snow_n, percentArea_snow_o, percentArea_snow_u, percent_snow_g, albedo_v_o, albedo_n_o, albedo_v_u, albedo_n_u, albedo_v_g, albedo_n_g, netRad_o, netRad_u, netRad_g, netRadLeaf_o_sunlit, netRadLeaf_o_shaded, netRadLeaf_u_sunlit, netRadLeaf_u_shaded, netShortRadLeaf_o_sunlit, netShortRadLeaf_o_shaded, netShortRadLeaf_u_sunlit, netShortRadLeaf_u_shaded)
end

function soilresp(Ccd, Cssd, Csmd, Cfsd, Cfmd, Csm, Cm, Cs, Cp, npp_yr, coef, soiltype, soilp, mid_res)
    ccall((:soilresp, libbeps), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cfloat, Ptr{Cdouble}, Cint, Ptr{Soil}, Ptr{Results}), Ccd, Cssd, Csmd, Cfsd, Cfmd, Csm, Cm, Cs, Cp, npp_yr, coef, soiltype, soilp, mid_res)
end

function readparam(lc, parameter1)
    ccall((:readparam, libbeps), Cvoid, (Cshort, Ptr{Cdouble}), lc, parameter1)
end

function lai2(stem_o, stem_u, LC, CosZs, lai_o, clumping, lai_u, lai_o_sunlit, lai_o_shaded, lai_u_sunlit, lai_u_shaded, PAI_o_sunlit, PAI_o_shaded, PAI_u_sunlit, PAI_u_shaded)
    ccall((:lai2, libbeps), Cvoid, (Cdouble, Cdouble, Cint, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), stem_o, stem_u, LC, CosZs, lai_o, clumping, lai_u, lai_o_sunlit, lai_o_shaded, lai_u_sunlit, lai_u_shaded, PAI_o_sunlit, PAI_o_shaded, PAI_u_sunlit, PAI_u_shaded)
end

function readcoef(lc, stxt, coef)
    ccall((:readcoef, libbeps), Cvoid, (Cshort, Cint, Ptr{Cdouble}), lc, stxt, coef)
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

function Leaf_Temperatures(Tair, slope, psychrometer, VPD_air, Cp_ca, Gw_o_sunlit, Gw_o_shaded, Gw_u_sunlit, Gw_u_shaded, Gww_o_sunlit, Gww_o_shaded, Gww_u_sunlit, Gww_u_shaded, Gh_o_sunlit, Gh_o_shaded, Gh_u_sunlit, Gh_u_shaded, Xcs_o, Xcl_o, Xcs_u, Xcl_u, radiation_o_sun, radiation_o_shaded, radiation_u_sun, radiation_u_shaded, Tc_o_sunlit, Tc_o_shaded, Tc_u_sunlit, Tc_u_shaded)
    ccall((:Leaf_Temperatures, libbeps), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), Tair, slope, psychrometer, VPD_air, Cp_ca, Gw_o_sunlit, Gw_o_shaded, Gw_u_sunlit, Gw_u_shaded, Gww_o_sunlit, Gww_o_shaded, Gww_u_sunlit, Gww_u_shaded, Gh_o_sunlit, Gh_o_shaded, Gh_u_sunlit, Gh_u_shaded, Xcs_o, Xcl_o, Xcs_u, Xcl_u, radiation_o_sun, radiation_o_shaded, radiation_u_sun, radiation_u_shaded, Tc_o_sunlit, Tc_o_shaded, Tc_u_sunlit, Tc_u_shaded)
end

function Leaf_Temperature(Tair, slope, psychrometer, VPD_air, Cp_ca, Gw, Gww, Gh, Xcs, Xcl, radiation)
    ccall((:Leaf_Temperature, libbeps), Cdouble, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble), Tair, slope, psychrometer, VPD_air, Cp_ca, Gw, Gww, Gh, Xcs, Xcl, radiation)
end

function sensible_heat(tempL_o_sunlit, tempL_o_shaded, tempL_u_sunlit, tempL_u_shaded, temp_g, temp_air, rh_air, Gheat_o_sunlit, Gheat_o_shaded, Gheat_u_sunlit, Gheat_u_shaded, Gheat_g, lai_o_sunlit, lai_o_shaded, lai_u_sunlit, lai_u_shaded, SH_o, SH_u, SH_g)
    ccall((:sensible_heat, libbeps), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), tempL_o_sunlit, tempL_o_shaded, tempL_u_sunlit, tempL_u_shaded, temp_g, temp_air, rh_air, Gheat_o_sunlit, Gheat_o_shaded, Gheat_u_sunlit, Gheat_u_shaded, Gheat_g, lai_o_sunlit, lai_o_shaded, lai_u_sunlit, lai_u_shaded, SH_o, SH_u, SH_g)
end

function transpiration(tempL_o_sunlit, tempL_o_shaded, tempL_u_sunlit, tempL_u_shaded, temp_air, rh_air, Gtrans_o_sunlit, Gtrans_o_shaded, Gtrans_u_sunlit, Gtrans_u_shaded, lai_o_sunlit, lai_o_shaded, lai_u_sunlit, lai_u_shaded, trans_o, trans_u)
    ccall((:transpiration, libbeps), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), tempL_o_sunlit, tempL_o_shaded, tempL_u_sunlit, tempL_u_shaded, temp_air, rh_air, Gtrans_o_sunlit, Gtrans_o_shaded, Gtrans_u_sunlit, Gtrans_u_shaded, lai_o_sunlit, lai_o_shaded, lai_u_sunlit, lai_u_shaded, trans_o, trans_u)
end

function evaporation_canopy(tempL_o_sunlit, tempL_o_shaded, tempL_u_sunlit, tempL_u_shaded, temp_air, rh_air, Gwater_o_sunlit, Gwater_o_shaded, Gwater_u_sunlit, Gwater_u_shaded, lai_o_sunlit, lai_o_shaded, lai_u_sunlit, lai_u_shaded, percent_water_o, percent_water_u, percent_snow_o, percent_snow_u, evapo_water_o, evapo_water_u, evapo_snow_o, evapo_snow_u)
    ccall((:evaporation_canopy, libbeps), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), tempL_o_sunlit, tempL_o_shaded, tempL_u_sunlit, tempL_u_shaded, temp_air, rh_air, Gwater_o_sunlit, Gwater_o_shaded, Gwater_u_sunlit, Gwater_u_shaded, lai_o_sunlit, lai_o_shaded, lai_u_sunlit, lai_u_shaded, percent_water_o, percent_water_u, percent_snow_o, percent_snow_u, evapo_water_o, evapo_water_u, evapo_snow_o, evapo_snow_u)
end

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
    ccall((:snowpack_stage1, libbeps), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), temp_air, precipitation, mass_snow_o_last, mass_snow_u_last, mass_snow_g_last, mass_snow_o, mass_snow_u, mass_snow_g, lai_o, lai_u, clumping, area_snow_o, area_snow_u, percent_snow_o, percent_snow_u, percent_snow_g, density_snow, depth_snow, albedo_v_snow, albedo_n_snow)
end

function snowpack_stage2(evapo_snow_o, evapo_snow_u, mass_snow_o, mass_snow_u)
    ccall((:snowpack_stage2, libbeps), Cvoid, (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), evapo_snow_o, evapo_snow_u, mass_snow_o, mass_snow_u)
end

function snowpack_stage3(temp_air, temp_snow, temp_snow_last, density_snow, depth_snow, depth_water, mass_snow_g)
    ccall((:snowpack_stage3, libbeps), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), temp_air, temp_snow, temp_snow_last, density_snow, depth_snow, depth_water, mass_snow_g)
end

const FW_VERSION = 1

const MAX_LAYERS = 10

const DEPTH_F = 6

const NOERROR = 0

const ERROR = 1

const PI = 3.1415926

const zero = 1.0e-10

const l_sta = 105

const l_end = 105

const p_sta = 101

const p_end = 101

const RTIMES = 24

const step = 3600

const kstep = 360

const kloop = 10

const layer = 5

const depth_f = 6

const CO2_air = 380

const rho_a = 1.292

# exports
const PREFIXES = ["CX", "clang_"]
for name in names(@__MODULE__; all=true), prefix in PREFIXES
    if startswith(string(name), prefix)
        @eval export $name
    end
end

end # module
