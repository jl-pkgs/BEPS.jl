# using BEPS
self(obj, args...) = obj
none(obj, args...) = nothing

export self, none

# abstract type AbstractHydroModel{FT} end
#<: AbstractHydroModel{FT}

# @with_kw_noshow
@with_kw mutable struct Soil{FT<:Real}
    # Properties belong to the whole soil profile
    flag::Int = 0# reserved for EnKF usage.
    n_layer::Int = 5 # the number of layers used in the model. Make sure n_layer <= MAX_LAYERS
    step_period::Int = 0

    MAX_LAYERS::Int = 10

    # Conditions on the top boundary
    Zp::FT = 0.0               # depth of ponded water on the ground surface, `[m]`
    Zsp::FT = 0.0              # snow depth, `[m]`
    r_rain_g::FT = 0.0         # the rainfall rate, un--on understory g--on ground surface, `[m/s]`
    soil_r::FT = 0.0           # soil surface resistance for water, discuss with Remi - an interface here
    r_drainage::FT = 0.0

    # Some variable used for soil
    r_root_decay::FT = 0.0    # decay_rate_of_root_distribution
    psi_min::FT = 0.0         # for fw
    alpha::FT = 0.0           # for fw
    f_soilwater::FT = 0.0

    # # /********************************************/
    # Properties belong to each soil horizon
    d_soil::Vector{FT} = zeros(FT, MAX_LAYERS)
    f_root::Vector{FT} = zeros(FT, MAX_LAYERS)       # root weight
    dt::Vector{FT} = zeros(FT, MAX_LAYERS)           # the weight calculated from soil_water_factor **re-calculate in the model

    # From read-param function
    thermal_cond::Vector{FT} = zeros(FT, MAX_LAYERS) # thermal conductivity. `[W m-1 K-1]`
    theta_vfc::Vector{FT} = zeros(FT, MAX_LAYERS)    # field capacity, `[m3 m-3]` (not used in this model. LHE. Feb. 01, 2013)
    theta_vwp::Vector{FT} = zeros(FT, MAX_LAYERS)    # wilt point    , `[m3 m-3]`
    fei::Vector{FT} = zeros(FT, MAX_LAYERS)          # porosity      , `[m3 m-3]`
    Ksat::Vector{FT} = zeros(FT, MAX_LAYERS)         # saturated hydraulic conductivity, `[m s-1]`
    psi_sat::Vector{FT} = zeros(FT, MAX_LAYERS)      # water potential at sat (`m`)
    b::Vector{FT} = zeros(FT, MAX_LAYERS)            # Cambell parameter b, `[-]`
    density_soil::Vector{FT} = zeros(FT, MAX_LAYERS) # soil bulk density  , `[kg m-3`]. LHE. Feb. 12, 2013.
    f_org::Vector{FT} = zeros(FT, MAX_LAYERS)        # volume fraction of organic matter, `[%]`.

    # Variables need to save
    ice_ratio::Vector{FT} = zeros(FT, MAX_LAYERS)    # the ratio of ice, (-)
    thetam::Vector{FT} = zeros(FT, MAX_LAYERS)
    thetam_prev::Vector{FT} = zeros(FT, MAX_LAYERS)  # soil water content, M3 m-3.

    # soil temperature in this layer, don't change because it is used in 
    # `soil_water_factor_v2`, and `UpdateSoil_Moisture`.
    # previous and current
    temp_soil_p::Vector{FT} = zeros(FT, MAX_LAYERS)
    temp_soil_c::Vector{FT} = zeros(FT, MAX_LAYERS)

    # Derived variables below:
    f_ice::Vector{FT} = zeros(FT, MAX_LAYERS)        # derived var.
    psib::Vector{FT} = zeros(FT, MAX_LAYERS)         # soil water suction at the bottom this layer
    psim::Vector{FT} = zeros(FT, MAX_LAYERS)         # soil water suction in this layer. 
    # Note: this variable can be derived from other parameters. LHE.

    thetab::Vector{FT} = zeros(FT, MAX_LAYERS)       # soil water content at the bottom of each layer
    r_waterflow::Vector{FT} = zeros(FT, MAX_LAYERS)  # the liquid water flow rates at the soil layer interfaces  
    # 'eg. 0,1,2..., represents the surface, the bottom of layer1, the bottom of layer2,...

    km::Vector{FT} = zeros(FT, MAX_LAYERS)
    Kb::Vector{FT} = zeros(FT, MAX_LAYERS) #the hydraulic conductivity of certain soil layer
    KK::Vector{FT} = zeros(FT, MAX_LAYERS)           # The average  conductivity of two soil layers.*/
    Cs::Vector{FT} = zeros(FT, MAX_LAYERS)

    # not used in gpp-only version. derived var.
    lambda::Vector{FT} = zeros(FT, MAX_LAYERS)       # thermal conductivity of each soil layer # /* ={0} by LHE */ 
    Ett::Vector{FT} = zeros(FT, MAX_LAYERS)          # ET in each layer. derived var

    # define a lambda_top for ice?
    G::Vector{FT} = zeros(FT, MAX_LAYERS)            # energy fluxes

end

function Soil(FT::DataType)
    Soil{FT}()
end

export Soil
