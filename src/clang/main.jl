var_o = zeros(41)
var_n = zeros(41)
v2last = zeros(41)
outp = zeros(3)
total = zeros(3)

p_soil = Soil()

LR = -200.0; #  -200.0 means no measured long-wave radiation, the value will be calculated later

parameter = readparam(landcover);      # n = 48
coef = readcoef(landcover, soil_type); # n = 48


for jday = 1:365
  # /***** re-calculate LAI & renew clump index *****/
  lai = lai_p[jday] * parameter[2] / clumping

  # Hour loop begin
  # for (rstep=0;rstep<24;rstep++)
  for rstep = 1:24
    flag = (jday == 1 && rstep == 1) ? 0 : 1

    # meteo->Srad = m_rad[jday-1][rstep];
    # meteo->temp = m_tem[jday-1][rstep];
    # meteo->rain = m_pre[jday-1][rstep];
    # meteo->wind = m_wind[jday-1][rstep];
    # meteo->LR = -200.0
    # tem = m_tem[jday-1][rstep];
    # hum = m_hum[jday-1][rstep];

    # Vapour pressure in mbar
    es = 0.46 * hum * (tem + 273.16) / 100
    esd = 6.1078 * exp((17.269 * tem) / (237.3 + tem))

    # Calculate relative humidity when reading in VPD
    rh = clamp(es / esd * 100, 0.0, 100.0) # relative humidity, %

    if (flag == 0)  # for 1st time step, to initialize var.
      # /***** initialize soil conditions, read soil parameters and set depth *****/
      Init_Soil_Parameters(landcover, soil_type, parameter[28], p_soil)
      p_soil.r_drainage = parameter[27]
      Init_Soil_Status(p_soil, st, tem, sw, snowdepth) # LHE

      # Initialize intermediate variables array
      var_o .= 0
      for i = 4:9
        var_o[i] = tem
      end
      for i = 10:15
        var_o[i] = p_soil.temp_soil_p[i-10]
      end
      for i = 22:27
        var_o[i] = p_soil.thetam_prev[i-22]
      end
      for i = 28:33
        var_o[i] = p_soil.ice_ratio[i-28]
      end
      # for (i=0;i<=40;i++)   var_o[i] = 0;
      # for (i=3;i<=8;i++)   var_o[i] = tem;
      # for(i=9;i<=14;i++) var_o[i] = p_soil->temp_soil_p[i-9];
      # for(i=21;i<=26;i++) var_o[i] = p_soil->thetam_prev[i-21];
      # for(i=27;i<=32;i++) var_o[i] = p_soil->ice_ratio[i-27];
    else   #  for other time steps, assigned from the temp variables array
      for i = 1:41
        var_o[i] = v2last[i]
      end
    end

    # /***** calculate cos_solar zenith angle Z *****/
    CosZs = s_coszs(jday, rstep, lat, lon)

    # /***** start simulation modules *****/
    #printf("%d, %d, %f\n", jday, rstep, p_soil->thetam_prev[0]); # in-process check
    inter_prg(jday, rstep, lai, clumping, parameter, meteo, CosZs, var_o, var_n, p_soil, mid_res)
    #printf("%d, %d, %f, %f\n", jday, rstep, p_soil->thetam_prev[0], p_soil->f_soilwater); # in-process check

    # Store updated variables array in temp array
    for i = 1:41
      v2last[i] = var_n[i]
    end
    # for (i=0;i<=40;i++)  

    # # if we don't care NEP, we can close this
    # # /***** plant respiration/NPP module *****/
    # temp_soil1=p_soil->temp_soil_c[1];
    # plantresp(landcover,mid_res,lai_yr,lai,tem,temp_soil1,CosZs);

    # # /***** soil respiration module *****/
    # soilresp(Ccd,Cssd,Csmd,Cfsd,Cfmd,Csm,Cm,Cs,Cp,nppyr,coef,soil_type,p_soil,mid_res);

    # /***** save data for output *****/
    # Hourly output
    outp[1] = mid_res.GPP
    outp[2] = mid_res.Trans + mid_res.Evap
    outp[3] = mid_res.NEP
    # outp[4]=mid_res->npp_o + mid_res->npp_u;
    # # Write hourly output to files
    # fprintf(outp_ptr,"%d %d gpp= %f tr= %f Ev= %f \n",jday,rstep,outp[1],outp[2],outp[3]);

    # # Sum of output
    total[1] = total[1] + outp[1]
    total[2] = total[2] + outp[2]
    total[3] = total[3] + outp[3]
  end # End of hourly loop
end # End of daily loop
