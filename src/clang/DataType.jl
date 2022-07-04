import Parameters: @with_kw, @with_kw_noshow

zero10() = tuple(zeros(10)...)

@with_kw mutable struct Soil
    flag::Cint = Cint(0)
    n_layer::Cint = Cint(5)
    step_period::Cint  = Cint(1)
    Zp::Cdouble = Cdouble(0)
    Zsp::Cdouble = Cdouble(0)
    r_rain_g::Cdouble = Cdouble(0)
    soil_r::Cdouble = Cdouble(0)
    r_drainage::Cdouble = Cdouble(0)
    r_root_decay::Cdouble = Cdouble(0)
    psi_min::Cdouble = Cdouble(0)
    alpha::Cdouble = Cdouble(0)
    f_soilwater::Cdouble = Cdouble(0)
    d_soil::NTuple{10, Cdouble} = zero10()
    f_root::NTuple{10, Cdouble} = zero10()
    dt::NTuple{10, Cdouble} = zero10()
    thermal_cond::NTuple{10, Cdouble} = zero10()
    theta_vfc::NTuple{10, Cdouble} = zero10()
    theta_vwp::NTuple{10, Cdouble} = zero10()
    fei::NTuple{10, Cdouble} = zero10()
    Ksat::NTuple{10, Cdouble} = zero10()
    psi_sat::NTuple{10, Cdouble} = zero10()
    b::NTuple{10, Cdouble} = zero10()
    density_soil::NTuple{10, Cdouble} = zero10()
    f_org::NTuple{10, Cdouble} = zero10()
    ice_ratio::NTuple{10, Cdouble} = zero10()
    thetam::NTuple{10, Cdouble} = zero10()
    thetam_prev::NTuple{10, Cdouble} = zero10()
    temp_soil_p::NTuple{10, Cdouble} = zero10()
    temp_soil_c::NTuple{10, Cdouble} = zero10()
    f_ice::NTuple{10, Cdouble} = zero10()
    psim::NTuple{10, Cdouble} = zero10()
    thetab::NTuple{10, Cdouble} = zero10()
    psib::NTuple{10, Cdouble} = zero10()
    r_waterflow::NTuple{10, Cdouble} = zero10()
    km::NTuple{10, Cdouble} = zero10()
    Kb::NTuple{10, Cdouble} = zero10()
    KK::NTuple{10, Cdouble} = zero10()
    Cs::NTuple{10, Cdouble} = zero10()
    lambda::NTuple{10, Cdouble} = zero10()
    Ett::NTuple{10, Cdouble} = zero10()
    G::NTuple{10, Cdouble} = zero10()
end


struct ClimateData
    Srad::Cdouble
    LR::Cdouble
    temp::Cdouble
    rh::Cdouble
    rain::Cdouble
    wind::Cdouble
    dr_o::Cdouble
    df_o::Cdouble
    dr_u::Cdouble
    df_u::Cdouble
end

struct Results
    gpp_o_sunlit::Cdouble
    gpp_u_sunlit::Cdouble
    gpp_o_shaded::Cdouble
    gpp_u_shaded::Cdouble
    plant_resp::Cdouble
    npp_o::Cdouble
    npp_u::Cdouble
    GPP::Cdouble
    NPP::Cdouble
    NEP::Cdouble
    soil_resp::Cdouble
    Net_Rad::Cdouble
    SH::Cdouble
    LH::Cdouble
    Trans::Cdouble
    Evap::Cdouble
end

struct Cpools
    Ccd::NTuple{3, Cdouble}
    Cssd::NTuple{3, Cdouble}
    Csmd::NTuple{3, Cdouble}
    Cfsd::NTuple{3, Cdouble}
    Cfmd::NTuple{3, Cdouble}
    Csm::NTuple{3, Cdouble}
    Cm::NTuple{3, Cdouble}
    Cs::NTuple{3, Cdouble}
    Cp::NTuple{3, Cdouble}
end
