export Init_Soil_Parameters, Init_Soil_Status, SoilRootFraction


"""
  Init_Soil_Parameters(landcover::Int, stxt::Int, r_root_decay::Float64, p::Soil)

# Example
```julia
x = Soil{Float64}(; MAX_LAYERS = 5)
Init_Soil_Parameters(1, 2, 0.1, x)
```
"""
function Init_Soil_Parameters(landcover::Int, stxt::Int, r_root_decay::Float64, p::Soil)
  p.n_layer = 5
  p.flag = 0
  p.step_period = 0
  p.soil_r = 0
  p.f_soilwater = 0

  if (landcover == 6 || landcover == 9)
    p.psi_min = 10.0
    p.alpha = 1.5
  else
    p.psi_min = 33.0
    p.alpha = 0.4
  end

  p.d_soil[1] = 0.05 ## /* depth_layer0 */
  p.d_soil[2] = 0.10 ## /* depth_layer1 */
  p.d_soil[3] = 0.20 ## /* depth_layer2 */
  p.d_soil[4] = 0.40 ## /* depth_layer3 */
  p.d_soil[5] = 1.25 ## /* depth_layer4 */

  p.r_root_decay = r_root_decay
  SoilRootFraction(p)

  p.density_soil[1] = 1300.0  # from flux tower.
  p.density_soil[2] = 1500.0
  p.density_soil[3] = 1517.0
  p.density_soil[4] = 1517.0
  p.density_soil[5] = 1517.0

  p.f_org[1] = 5
  p.f_org[2] = 2
  p.f_org[3] = 1
  p.f_org[4] = 1
  p.f_org[5] = 0.3

  if stxt == 1
    #sand
    p.b[1] = 1.7
    p.b[2] = 1.9
    p.b[3] = 2.1
    p.b[4] = 2.3
    p.b[5] = 2.5
    p.Ksat[1] = 0.000058
    p.Ksat[2] = 0.000052
    p.Ksat[3] = 0.000046
    p.Ksat[4] = 0.000035
    p.Ksat[5] = 0.000010 # /* saturated hydraulic conductivity */
    p.fei[1] = 0.437
    p.fei[2] = 0.437
    p.fei[3] = 0.437
    p.fei[4] = 0.437
    p.fei[5] = 0.437 # /* porosity */
    p.theta_vfc[1] = 0.09
    p.theta_vfc[2] = 0.09
    p.theta_vfc[3] = 0.09
    p.theta_vfc[4] = 0.09
    p.theta_vfc[5] = 0.09 # /* field capacity */
    p.theta_vwp[1] = 0.03
    p.theta_vwp[2] = 0.03
    p.theta_vwp[3] = 0.03
    p.theta_vwp[4] = 0.03
    p.theta_vwp[5] = 0.03 # /* wilt point*/
    p.thermal_cond[1] = 8.6
    p.thermal_cond[2] = 8.6
    p.thermal_cond[3] = 8.6
    p.thermal_cond[4] = 8.6
    p.thermal_cond[5] = 8.6 # /* thermal conductivity */
    p.psi_sat[1] = 0.07
    p.psi_sat[2] = 0.08
    p.psi_sat[3] = 0.09
    p.psi_sat[4] = 0.10
    p.psi_sat[5] = 0.12 # /* water potential at sat */
  # break;
  elseif stxt == 2
    # loamy sand
    p.b[1] = 2.1
    p.b[2] = 2.3
    p.b[3] = 2.5
    p.b[4] = 2.7
    p.b[5] = 2.9
    p.Ksat[1] = 0.000017
    p.Ksat[2] = 0.000015
    p.Ksat[3] = 0.000014
    p.Ksat[4] = 0.000010
    p.Ksat[5] = 0.000003
    p.fei[1] = 0.437
    p.fei[2] = 0.437
    p.fei[3] = 0.437
    p.fei[4] = 0.437
    p.fei[5] = 0.437
    p.theta_vfc[1] = 0.21
    p.theta_vfc[2] = 0.21
    p.theta_vfc[3] = 0.21
    p.theta_vfc[4] = 0.21
    p.theta_vfc[5] = 0.21
    p.theta_vwp[1] = 0.06
    p.theta_vwp[2] = 0.06
    p.theta_vwp[3] = 0.06
    p.theta_vwp[4] = 0.06
    p.theta_vwp[5] = 0.06
    p.thermal_cond[1] = 8.3
    p.thermal_cond[2] = 8.3
    p.thermal_cond[3] = 8.3
    p.thermal_cond[4] = 8.3
    p.thermal_cond[5] = 8.3
    p.psi_sat[1] = 0.09
    p.psi_sat[2] = 0.10
    p.psi_sat[3] = 0.11
    p.psi_sat[4] = 0.12
    p.psi_sat[5] = 0.14
    # break;
  elseif stxt == 3
    # case 3:  # sandy loam
    p.b[1] = 3.1
    p.b[2] = 3.3
    p.b[3] = 3.5
    p.b[4] = 3.7
    p.b[5] = 3.9
    p.Ksat[1] = 0.0000072
    p.Ksat[2] = 0.00000648
    p.Ksat[3] = 0.00000576
    p.Ksat[4] = 0.00000432
    p.Ksat[5] = 0.00000144
    p.fei[1] = 0.453
    p.fei[2] = 0.453
    p.fei[3] = 0.453
    p.fei[4] = 0.453
    p.fei[5] = 0.453
    p.theta_vfc[1] = 0.21
    p.theta_vfc[2] = 0.21
    p.theta_vfc[3] = 0.21
    p.theta_vfc[4] = 0.21
    p.theta_vfc[5] = 0.21
    p.theta_vwp[1] = 0.10
    p.theta_vwp[2] = 0.10
    p.theta_vwp[3] = 0.10
    p.theta_vwp[4] = 0.10
    p.theta_vwp[5] = 0.10
    p.thermal_cond[1] = 8.0
    p.thermal_cond[2] = 8.0
    p.thermal_cond[3] = 8.0
    p.thermal_cond[4] = 8.0
    p.thermal_cond[5] = 8.0
    p.psi_sat[1] = 0.15
    p.psi_sat[2] = 0.16
    p.psi_sat[3] = 0.17
    p.psi_sat[4] = 0.18
    p.psi_sat[5] = 0.20
    # break;
  elseif stxt == 4
    # case 4:  # loam
    p.b[1] = 4.5
    p.b[2] = 4.7
    p.b[3] = 4.9
    p.b[4] = 5.1
    p.b[5] = 5.3
    p.Ksat[1] = 0.0000037
    p.Ksat[2] = 0.0000033
    p.Ksat[3] = 0.00000296
    p.Ksat[4] = 0.00000222
    p.Ksat[5] = 0.00000074
    p.fei[1] = 0.463
    p.fei[2] = 0.463
    p.fei[3] = 0.463
    p.fei[4] = 0.463
    p.fei[5] = 0.463
    p.theta_vfc[1] = 0.27
    p.theta_vfc[2] = 0.27
    p.theta_vfc[3] = 0.27
    p.theta_vfc[4] = 0.27
    p.theta_vfc[5] = 0.27
    p.theta_vwp[1] = 0.12
    p.theta_vwp[2] = 0.12
    p.theta_vwp[3] = 0.12
    p.theta_vwp[4] = 0.12
    p.theta_vwp[5] = 0.12
    p.thermal_cond[1] = 7.0
    p.thermal_cond[2] = 7.0
    p.thermal_cond[3] = 7.0
    p.thermal_cond[4] = 7.0
    p.thermal_cond[5] = 7.0
    p.psi_sat[1] = 0.11
    p.psi_sat[2] = 0.12
    p.psi_sat[3] = 0.13
    p.psi_sat[4] = 0.14
    p.psi_sat[5] = 0.16
    # break;
  elseif stxt == 5
    # case 5:  # silty loam
    p.b[1] = 4.7
    p.b[2] = 4.9
    p.b[3] = 5.1
    p.b[4] = 5.3
    p.b[5] = 5.5
    p.Ksat[1] = 0.0000019
    p.Ksat[2] = 0.0000017
    p.Ksat[3] = 0.00000152
    p.Ksat[4] = 0.00000114
    p.Ksat[5] = 0.00000038
    p.fei[1] = 0.501
    p.fei[2] = 0.501
    p.fei[3] = 0.501
    p.fei[4] = 0.501
    p.fei[5] = 0.501
    p.theta_vfc[1] = 0.33
    p.theta_vfc[2] = 0.33
    p.theta_vfc[3] = 0.33
    p.theta_vfc[4] = 0.33
    p.theta_vfc[5] = 0.33
    p.theta_vwp[1] = 0.13
    p.theta_vwp[2] = 0.13
    p.theta_vwp[3] = 0.13
    p.theta_vwp[4] = 0.13
    p.theta_vwp[5] = 0.13
    p.thermal_cond[1] = 6.3
    p.thermal_cond[2] = 6.3
    p.thermal_cond[3] = 6.3
    p.thermal_cond[4] = 6.3
    p.thermal_cond[5] = 6.3
    p.psi_sat[1] = 0.21
    p.psi_sat[2] = 0.22
    p.psi_sat[3] = 0.23
    p.psi_sat[4] = 0.24
    p.psi_sat[5] = 0.26
    # break;
  elseif stxt == 6
    # case 6:  # sandy clay loam
    p.b[1] = 4.0
    p.b[2] = 4.2
    p.b[3] = 4.4
    p.b[4] = 4.6
    p.b[5] = 4.8
    p.Ksat[1] = 0.0000012
    p.Ksat[2] = 0.00000108
    p.Ksat[3] = 0.0000096
    p.Ksat[4] = 0.0000072
    p.Ksat[5] = 0.0000024
    p.fei[1] = 0.398
    p.fei[2] = 0.398
    p.fei[3] = 0.398
    p.fei[4] = 0.398
    p.fei[5] = 0.398
    p.theta_vfc[1] = 0.26
    p.theta_vfc[2] = 0.26
    p.theta_vfc[3] = 0.26
    p.theta_vfc[4] = 0.26
    p.theta_vfc[5] = 0.26
    p.theta_vwp[1] = 0.15
    p.theta_vwp[2] = 0.15
    p.theta_vwp[3] = 0.15
    p.theta_vwp[4] = 0.15
    p.theta_vwp[5] = 0.15
    p.thermal_cond[1] = 7.0
    p.thermal_cond[2] = 7.0
    p.thermal_cond[3] = 7.0
    p.thermal_cond[4] = 7.0
    p.thermal_cond[5] = 7.0
    p.psi_sat[1] = 0.28
    p.psi_sat[2] = 0.29
    p.psi_sat[3] = 0.30
    p.psi_sat[4] = 0.31
    p.psi_sat[5] = 0.33
    # break;
  elseif stxt == 7
    # case 7:  # clay loam
    p.b[1] = 5.2
    p.b[2] = 5.4
    p.b[3] = 5.6
    p.b[4] = 5.8
    p.b[5] = 6.0
    p.Ksat[1] = 0.00000064
    p.Ksat[2] = 0.00000058
    p.Ksat[3] = 0.00000051
    p.Ksat[4] = 0.00000038
    p.Ksat[5] = 0.00000013
    p.fei[1] = 0.464
    p.fei[2] = 0.464
    p.fei[3] = 0.464
    p.fei[4] = 0.464
    p.fei[5] = 0.464
    p.theta_vfc[1] = 0.32
    p.theta_vfc[2] = 0.32
    p.theta_vfc[3] = 0.32
    p.theta_vfc[4] = 0.32
    p.theta_vfc[5] = 0.32
    p.theta_vwp[1] = 0.20
    p.theta_vwp[2] = 0.20
    p.theta_vwp[3] = 0.20
    p.theta_vwp[4] = 0.20
    p.theta_vwp[5] = 0.20
    p.thermal_cond[1] = 5.8
    p.thermal_cond[2] = 5.8
    p.thermal_cond[3] = 5.7
    p.thermal_cond[4] = 5.8
    p.thermal_cond[5] = 5.8
    p.psi_sat[1] = 0.26
    p.psi_sat[2] = 0.27
    p.psi_sat[3] = 0.28
    p.psi_sat[4] = 0.29
    p.psi_sat[5] = 0.31
    # break;
  elseif stxt == 8
    # case 8:  # silty clay loam
    p.b[1] = 6.6
    p.b[2] = 6.8
    p.b[3] = 7.0
    p.b[4] = 7.2
    p.b[5] = 7.4
    p.Ksat[1] = 0.00000042
    p.Ksat[2] = 0.00000038
    p.Ksat[3] = 0.00000034
    p.Ksat[4] = 0.000000252
    p.Ksat[5] = 0.000000084
    p.fei[1] = 0.471
    p.fei[2] = 0.471
    p.fei[3] = 0.471
    p.fei[4] = 0.471
    p.fei[5] = 0.471
    p.theta_vfc[1] = 0.37
    p.theta_vfc[2] = 0.37
    p.theta_vfc[3] = 0.37
    p.theta_vfc[4] = 0.37
    p.theta_vfc[5] = 0.37
    p.theta_vwp[1] = 0.32
    p.theta_vwp[2] = 0.32
    p.theta_vwp[3] = 0.32
    p.theta_vwp[4] = 0.32
    p.theta_vwp[5] = 0.32
    p.thermal_cond[1] = 4.2
    p.thermal_cond[2] = 4.2
    p.thermal_cond[3] = 4.2
    p.thermal_cond[4] = 4.2
    p.thermal_cond[5] = 4.2
    p.psi_sat[1] = 0.33
    p.psi_sat[2] = 0.34
    p.psi_sat[3] = 0.35
    p.psi_sat[4] = 0.36
    p.psi_sat[5] = 0.38
    # break;

  elseif stxt == 9
    # case 9:  # sandy clay
    p.b[1] = 6.0
    p.b[2] = 6.2
    p.b[3] = 6.4
    p.b[4] = 6.6
    p.b[5] = 6.8
    p.Ksat[1] = 0.00000033
    p.Ksat[2] = 0.0000003
    p.Ksat[3] = 0.000000264
    p.Ksat[4] = 0.000000198
    p.Ksat[5] = 0.000000066
    p.fei[1] = 0.430
    p.fei[2] = 0.430
    p.fei[3] = 0.430
    p.fei[4] = 0.430
    p.fei[5] = 0.430
    p.theta_vfc[1] = 0.34
    p.theta_vfc[2] = 0.34
    p.theta_vfc[3] = 0.34
    p.theta_vfc[4] = 0.34
    p.theta_vfc[5] = 0.34
    p.theta_vwp[1] = 0.24
    p.theta_vwp[2] = 0.24
    p.theta_vwp[3] = 0.24
    p.theta_vwp[4] = 0.24
    p.theta_vwp[5] = 0.24
    p.thermal_cond[1] = 6.3
    p.thermal_cond[2] = 6.3
    p.thermal_cond[3] = 6.3
    p.thermal_cond[4] = 6.3
    p.thermal_cond[5] = 6.3
    p.psi_sat[1] = 0.29
    p.psi_sat[2] = 0.30
    p.psi_sat[3] = 0.31
    p.psi_sat[4] = 0.32
    p.psi_sat[5] = 0.34
    # break;
  elseif stxt == 10
    # case 10:  # silty clay
    p.b[1] = 7.9
    p.b[2] = 8.1
    p.b[3] = 8.3
    p.b[4] = 8.5
    p.b[5] = 8.7
    p.Ksat[1] = 0.00000025
    p.Ksat[2] = 0.000000225
    p.Ksat[3] = 0.0000002
    p.Ksat[4] = 0.00000015
    p.Ksat[5] = 0.00000005
    p.fei[1] = 0.479
    p.fei[2] = 0.479
    p.fei[3] = 0.479
    p.fei[4] = 0.479
    p.fei[5] = 0.479
    p.theta_vfc[1] = 0.39
    p.theta_vfc[2] = 0.39
    p.theta_vfc[3] = 0.39
    p.theta_vfc[4] = 0.39
    p.theta_vfc[5] = 0.39
    p.theta_vwp[1] = 0.25
    p.theta_vwp[2] = 0.25
    p.theta_vwp[3] = 0.25
    p.theta_vwp[4] = 0.25
    p.theta_vwp[5] = 0.25
    p.thermal_cond[1] = 4.0
    p.thermal_cond[2] = 4.0
    p.thermal_cond[3] = 4.0
    p.thermal_cond[4] = 4.0
    p.thermal_cond[5] = 4.0
    p.psi_sat[1] = 0.34
    p.psi_sat[2] = 0.35
    p.psi_sat[3] = 0.36
    p.psi_sat[4] = 0.37
    p.psi_sat[5] = 0.39
    # break;
  elseif stxt == 11
    # case 11:  # clay
    p.b[1] = 7.6
    p.b[2] = 7.8
    p.b[3] = 8.0
    p.b[4] = 8.2
    p.b[5] = 8.4
    p.Ksat[1] = 0.00000017
    p.Ksat[2] = 0.000000153
    p.Ksat[3] = 0.000000136
    p.Ksat[4] = 0.000000102
    p.Ksat[5] = 0.000000034
    p.fei[1] = 0.475
    p.fei[2] = 0.475
    p.fei[3] = 0.475
    p.fei[4] = 0.475
    p.fei[5] = 0.475
    p.theta_vfc[1] = 0.40
    p.theta_vfc[2] = 0.40
    p.theta_vfc[3] = 0.40
    p.theta_vfc[4] = 0.40
    p.theta_vfc[5] = 0.40
    p.theta_vwp[1] = 0.27
    p.theta_vwp[2] = 0.27
    p.theta_vwp[3] = 0.27
    p.theta_vwp[4] = 0.27
    p.theta_vwp[5] = 0.27
    p.thermal_cond[1] = 4.4
    p.thermal_cond[2] = 4.4
    p.thermal_cond[3] = 4.4
    p.thermal_cond[4] = 4.4
    p.thermal_cond[5] = 4.4
    p.psi_sat[1] = 0.37
    p.psi_sat[2] = 0.38
    p.psi_sat[3] = 0.39
    p.psi_sat[4] = 0.40
    p.psi_sat[5] = 0.42
    # break;
  else
    # default:
    p.b[1] = 7.6
    p.b[2] = 7.8
    p.b[3] = 8.0
    p.b[4] = 8.2
    p.b[5] = 8.4
    p.Ksat[1] = 0.00000017
    p.Ksat[2] = 0.000000153
    p.Ksat[3] = 0.000000136
    p.Ksat[4] = 0.000000102
    p.Ksat[5] = 0.000000034
    p.fei[1] = 0.475
    p.fei[2] = 0.475
    p.fei[3] = 0.475
    p.fei[4] = 0.475
    p.fei[5] = 0.475
    p.theta_vfc[1] = 0.40
    p.theta_vfc[2] = 0.40
    p.theta_vfc[3] = 0.40
    p.theta_vfc[4] = 0.40
    p.theta_vfc[5] = 0.40
    p.theta_vwp[1] = 0.27
    p.theta_vwp[2] = 0.27
    p.theta_vwp[3] = 0.27
    p.theta_vwp[4] = 0.27
    p.theta_vwp[5] = 0.27
    p.thermal_cond[1] = 4.4
    p.thermal_cond[2] = 4.4
    p.thermal_cond[3] = 4.4
    p.thermal_cond[4] = 4.4
    p.thermal_cond[5] = 4.4
    p.psi_sat[1] = 0.37
    p.psi_sat[2] = 0.38
    p.psi_sat[3] = 0.39
    p.psi_sat[4] = 0.40
    p.psi_sat[5] = 0.42
  end
end

# #/ @brief Function to initialize the soil status:
# #/        soil temperature and moisture for each layer,
# #/        ponded water, snow depth, et al.
# #/ @param  p          Soil struct variable
# #/ @param  Tsoil      soil temperature
# #/ @param  Tair       air temperature
# #/ @param  Ms         soil water content
# #/ @param  snowdepth  snow depth
function Init_Soil_Status(p::Soil, Tsoil::FT, Tair::FT, Ms::FT, snowdepth::FT) where {FT<:Real}
  d_t = Tsoil - Tair

  p.Zp = 0.0        # /* depth of ponded water on the surface*/
  p.Zsp = snowdepth # /*Snow depth;*/
  p.r_rain_g = 0.0  # /*the rainfall rate, un--on understory g--on ground surface  m/s */

  d_t = clamp(d_t, -5.0, 5.0)
  # if (d_t > 5.0); d_t = 5.0; end
  # if (d_t < -5.0); d_t = -5.0; end

  p.temp_soil_c[1] = Tair + 0.4 * d_t
  p.temp_soil_c[2] = Tair + 0.5 * d_t
  p.temp_soil_c[3] = Tair + d_t
  p.temp_soil_c[4] = Tair + 1.2 * d_t
  p.temp_soil_c[5] = Tair + 1.4 * d_t

  p.temp_soil_p[1] = Tair + 0.4 * d_t
  p.temp_soil_p[2] = Tair + 0.5 * d_t
  p.temp_soil_p[3] = Tair + d_t
  p.temp_soil_p[4] = Tair + 1.2 * d_t
  p.temp_soil_p[5] = Tair + 1.4 * d_t

  #p.thetam[0] = 0.75*Ms;
  p.thetam[1] = 0.8 * Ms
  p.thetam[2] = Ms
  p.thetam[3] = 1.05 * Ms
  p.thetam[4] = 1.10 * Ms
  p.thetam[5] = 1.15 * Ms

  #p.thetam_prev[0] = 0.75*Ms;
  p.thetam_prev[1] = 0.8 * Ms
  p.thetam_prev[2] = Ms
  p.thetam_prev[3] = 1.05 * Ms
  p.thetam_prev[4] = 1.10 * Ms
  p.thetam_prev[5] = 1.15 * Ms

  # for (i = 0; i < p.n_layer; i++) {
  for i = 1:p.n_layer
    if (p.temp_soil_c[i] < -1.0)
      p.ice_ratio[i] = 1.0
    elseif (p.temp_soil_c[i] > 0)
      p.ice_ratio[i] = 0
    else
      p.ice_ratio[i] = (0 - p.temp_soil_c[i]) / 1.0
    end
  end
end

function SoilRootFraction(soil::Soil)
  # double cum_depth[soil.MAX_LAYERS];
  cum_depth = zeros(soil.MAX_LAYERS)
  n = soil.n_layer
  # for the 0 layer
  cum_depth[1] = soil.d_soil[1]
  soil.f_root[1] = 1 - pow(soil.r_root_decay, cum_depth[1] * 100)

  for i = 2:n-1
    cum_depth[i] = cum_depth[i-1] + soil.d_soil[i]
    soil.f_root[i] = pow(soil.r_root_decay, cum_depth[i-1] * 100) - pow(soil.r_root_decay, cum_depth[i] * 100)
  end
  # for the last layer. Put all reminding roots to the last layer.
  soil.f_root[n] = pow(soil.r_root_decay, cum_depth[n-1] * 100)
end
