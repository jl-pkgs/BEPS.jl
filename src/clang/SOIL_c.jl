include("DataType.jl")


function SoilRootFraction(soil)
    ccall((:SoilRootFraction, libbeps), Cvoid, (Ptr{Soil},), soil)
end

function Init_Soil_Parameters(landcover::Int, stxt::Int, r_root_decay::Real, p::Soil)
    ccall((:Init_Soil_Parameters, libbeps), Cvoid, (Cint, Cint, Cdouble, Ptr{Soil}), 
    Cint(landcover), Cint(stxt), Cfloat(r_root_decay), Ref(p))
end

function Init_Soil_Status(p, Tsoil, Tair, Ms, snowdepth)
    ccall((:Init_Soil_Status, libbeps), Cvoid, (Ptr{Soil}, Cdouble, Cdouble, Cdouble, Cdouble), p, Tsoil, Tair, Ms, snowdepth)
end

function soil_water_factor_v2(p)
    ccall((:soil_water_factor_v2, libbeps), Cvoid, (Ptr{Soil},), p)
end

function Soil_Water_Uptake(p, Trans_o, Trans_u, Evap_soil)
    ccall((:Soil_Water_Uptake, libbeps), Cvoid, (Ptr{Soil}, Cdouble, Cdouble, Cdouble), p, Trans_o, Trans_u, Evap_soil)
end

function UpdateSoilLambda(soil)
    ccall((:UpdateSoilLambda, libbeps), Cvoid, (Ptr{Soil},), soil)
end

function init_soil_parameter(T_USDA, S_USDA, Ref_Depth, T_Density, S_Density, T_OC, S_OC, soil)
    ccall((:init_soil_parameter, libbeps), Cvoid, (Cuchar, Cuchar, Cuchar, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Soil}), T_USDA, S_USDA, Ref_Depth, T_Density, S_Density, T_OC, S_OC, soil)
end

function Update_Cs(p)
    ccall((:Update_Cs, libbeps), Cvoid, (Ptr{Soil},), p)
end

function Update_ice_ratio(p)
    ccall((:Update_ice_ratio, libbeps), Cvoid, (Ptr{Soil},), p)
end

function UpdateSoilThermalConductivity(p)
    ccall((:UpdateSoilThermalConductivity, libbeps), Cvoid, (Ptr{Soil},), p)
end

function UpdateHeatFlux(p, Xg_snow, lambda_snow, Tsn0, Tair_annual_mean, peroid_in_seconds)
    ccall((:UpdateHeatFlux, libbeps), Cvoid, (Ptr{Soil}, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble), p, Xg_snow, lambda_snow, Tsn0, Tair_annual_mean, peroid_in_seconds)
end

function UpdateSoilMoisture(p, peroid_in_seconds)
    ccall((:UpdateSoilMoisture, libbeps), Cvoid, (Ptr{Soil}, Cdouble), p, peroid_in_seconds)
end
