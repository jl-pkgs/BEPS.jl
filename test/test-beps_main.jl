using BEPS
using BEPS.beps_c
# using UnPack
# using DelimitedFiles: readdlm
# include("main_pkgs.jl")

begin
  f = "examples/input/p1_meteo.txt"
  d = fread(f)
  lai = readdlm("examples/input/p1_lai.txt")[:]

  par = (lon=120.5, lat=30.5, landcover=25, clumping=0.85,
    soil_type=8, Tsoil=2.2, soilwater=0.4115, snowdepth=0)
  # @unpack lon, lat, landcover, clumping, soil_type, Tsoil, soilwater, snowdepth = par

  # 2.1  990  3.4124  0.6607  0.0177  0.97491  0.0213  0.0531  0.9715  12.5484  8.8816      
  # LAI_yr, ann_NPP,  ccd,  cssd,  csmd,  cfsd,  cfmd,  csm,  cm,  cs,  cp

  # 120.5  30.5  25  0.85  8  2.2  0.4115  0.0      
  # long., lat., LC, CI, soiltxt, Tsoil, soilwater,snow-dp
  ""
end

@time df_out = besp_main(d, lai, par);
sum(df_out)

r1 = r2 = Ref(1.0)

# :gpp_o_sunlit => 41704.09054753295
# :gpp_u_sunlit => 285.1932406544977
# :gpp_o_shaded => 12720.54113242912
# :gpp_u_shaded => 135.13475998737647
# :plant_resp => 0.0
# :npp_o => 0.0
# :npp_u => 0.0
# :GPP => 2369.3022582020903
# :NPP => 0.0
# :NEP => 0.0
# :soil_resp => 0.0
# :Net_Rad => 1.299620168997007e6
# :SH => 232161.16420363513
# :LH => 811896.0668881282
# :Trans => 444.02437414775443
# :Evap => 748.8880942489606

# fwrite(df_out, "p1_out.csv")
