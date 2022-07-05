include("DataType.jl")


function SoilRootFraction(p::Soil)
    ccall((:SoilRootFraction, libbeps), Cvoid, (Ptr{Soil},), Ref(p))
end

function Init_Soil_Parameters(landcover::Int, stxt::Int, r_root_decay::Real, p::Soil)
    ccall((:Init_Soil_Parameters, libbeps), Cvoid,
        (Cint, Cint, Cdouble, Ptr{Soil}),
        Cint(landcover), Cint(stxt), Cfloat(r_root_decay), Ref(p))
end

function Init_Soil_Status(p::Soil, Tsoil, Tair, Ms, snowdepth)
    ccall((:Init_Soil_Status, libbeps), Cvoid,
        (Ptr{Soil}, Cdouble, Cdouble, Cdouble, Cdouble),
        Ref(p), Tsoil, Tair, Ms, snowdepth)
end

function soil_water_factor_v2(p::Soil)
    ccall((:soil_water_factor_v2, libbeps), Cvoid, (Ptr{Soil},), Ref(p))
end

function Soil_Water_Uptake(p::Soil, Trans_o, Trans_u, Evap_soil)
    ccall((:Soil_Water_Uptake, libbeps), Cvoid,
        (Ptr{Soil}, Cdouble, Cdouble, Cdouble),
        Ref(p), Trans_o, Trans_u, Evap_soil)
end

function UpdateSoilLambda(p::Soil)
    ccall((:UpdateSoilLambda, libbeps), Cvoid, (Ptr{Soil},), Ref(p))
end

function init_soil_parameter(T_USDA, S_USDA, Ref_Depth, T_Density, S_Density, T_OC, S_OC, p::Soil)
    ccall((:init_soil_parameter, libbeps), Cvoid,
        (Cuchar, Cuchar, Cuchar, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Soil}),
        T_USDA, S_USDA, Ref_Depth, T_Density, S_Density, T_OC, S_OC, Ref(p))
end

function Update_Cs(p::Soil)
    ccall((:Update_Cs, libbeps), Cvoid, (Ptr{Soil},), Ref(p))
end

function Update_ice_ratio(p::Soil)
    ccall((:Update_ice_ratio, libbeps), Cvoid, (Ptr{Soil},), Ref(p))
end

function UpdateSoilThermalConductivity(p::Soil)
    ccall((:UpdateSoilThermalConductivity, libbeps), Cvoid, (Ptr{Soil},), Ref(p))
end

function UpdateHeatFlux(p::Soil, Xg_snow, lambda_snow, Tsn0, Tair_annual_mean, peroid_in_seconds)
    ccall((:UpdateHeatFlux, libbeps), Cvoid,
        (Ptr{Soil}, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble),
        Ref(p), Xg_snow, lambda_snow, Tsn0, Tair_annual_mean, peroid_in_seconds)
end

function UpdateSoilMoisture(p::Soil, peroid_in_seconds)
    ccall((:UpdateSoilMoisture, libbeps), Cvoid, (Ptr{Soil}, Cdouble),
        Ref(p), peroid_in_seconds)
end
