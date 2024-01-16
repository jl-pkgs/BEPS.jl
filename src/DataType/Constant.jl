const FW_VERSION = 1

const MAX_LAYERS = 10

const DEPTH_F = 6

const NOERROR = 0

const ERROR = 1

const PI::Float64 = 3.1415926

const zero::Float64 = 1.0e-10

const l_sta = 105

const l_end = 105

const p_sta = 101

const p_end = 101

const RTIMES = 24

const step::Float64 = 3600.0

const kstep::Float64 = 360.0

const kloop = 10

const MAX_Loop = 11

const layer = 5

const depth_f = 6

const CO2_air::Float64 = 380.0    # ppm

const rho_a::Float64 = 1.292    # density of air, kg m-3

const rho_w::Float64 = 1025.0;  # density of water, kg m-3

const Cpd::Float64 = 1004.65;

const Lv_solid::Float64 = 2.83 * 1000000;  # the latent heat of evaporation from solid (snow/ice) at air temperature=Ta, in j+kkk/kg

## DB ----------------

const pi180::Float64 = 0.017453292  # pi divided by 180, radians per degree
const pi9::Float64 = 2.864788976
const PI2::Float64 = 6.283185307  # 2 times pi
const PI4::Float64 = 12.5663706

###
const rugc::Float64 = 8.314  # J mole-1 K-1
# const vcopt = 73.0   # carboxylation rate at optimal temperature, umol m-2 s-1
# const jmopt = 170.0  # electron transport rate at optimal temperature, umol m-2 s-1
const rd25::Float64 = 0.34  # dark respiration at 25 C, rd25= 0.34 umol m-2 s-1

# Universal gas constant
const rgc1000::Float64 = 8314.0  # gas constant times 1000.

# Consts for Photosynthesis model and kinetic equations.
# for Vcmax and Jmax.  Taken from Harley and Baldocchi (1995, PCE)
# const hkin::Float64 = 200000.0  # enthalpy term, J mol-1
const skin::Float64 = 710.0     # entropy term, J K-1 mol-1
const ejm::Float64 = 55000.0    # activation energy for electron transport, J mol-1
const evc::Float64 = 55000.0    # activation energy for carboxylation, J mol-1

# Enzyme constants & partial pressure of O2 and CO2
# Michaelis-Menten K values. From survey of literature.

const kc25::Float64 = 274.6  # kinetic coef for CO2 at 25 C, microbars
const ko25::Float64 = 419.8  # kinetic coef for O2 at 25C,  millibars
const o2::Float64 = 210.0    # oxygen concentration  mmol mol-1

# tau is computed on the basis of the Specificity factor (102.33)
# times Kco2/Kh2o (28.38) to convert for value in solution
# to that based in air/
# The old value was 2321.1.

# New value for Quercus robor from Balaguer et al. 1996
# Similar number from Dreyer et al. 2001, Tree Physiol, tau= 2710

const tau25::Float64 = 2904.12   #  tau coefficient
#  Arrhenius constants
#  Eact for Michaelis-Menten const. for KC, KO and dark respiration
#  These values are from Harley
const ekc::Float64 = 80500.0     # Activation energy for K of CO2; J mol-1
const eko::Float64 = 14500.0     # Activation energy for K of O2, J mol-1
const erd::Float64 = 38000.0     # activation energy for dark respiration, eg Q10=2
const ektau::Float64 = -29000.0  # J mol-1 (Jordan and Ogren, 1984)
const tk_25::Float64 = 298.16    # absolute temperature at 25 C
const toptvc::Float64 = 301.0    # optimum temperature for maximum carboxylation
const toptjm::Float64 = 301.0    # optimum temperature for maximum electron transport
const eabole::Float64 = 45162.0    # activation energy for bole respiration for Q10 = 2.02
