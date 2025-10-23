_types = ["ENF", "DNF", "DBF", "EBF", "Shrub-SH", "C4", "default"]
_codes = [1, 2, 6, 9, 13, 40, -1]

"""
    readVegParam(lc::Int=1)
# Arguments
- `lc`: [1, 2, 6, 9, 13, 40, -1]
"""
function readVegParam(lc::Int=1)
  param = zeros(Float64, 48)
  param[5] = lc

  if lc == 1
    # Conifer Evergreen Forest (ENF) parameters
    param[3] = 0.62                     # clumping_index  
    param[9] = 4.5                      # LAI_max_overstory                     3.3
    param[10] = 2.4                     # LAI_max_understory                    2.4
    param[17] = 1.33                    # z00 roughness length, for heat = 0.1 canopy height = 1.5 
    param[19] = 35                      # mass_overstory(kg/m^2)                25
    param[21] = 10                      # mass_understory(kg/m^2)               10
    param[22] = 0.6                     # root_depth(m)                         0.8
    param[23] = 0.035                   # albedo_canopy_vis                     0.09
    param[24] = 0.23                    # albedo_c_nir_observed                 0.19
    param[28] = 0.95                    # decay_rate_of_root_distribution       0.97
    param[29] = 150                     # minimum_stomatal_resistance(s/m)      100
    param[30] = 20                      # canopy_height                         23
    param[31] = 3                       # height of understory                  3

    # the following parameters were added by MEM for the Ball and Berry calculations
    param[34] = 8                       # k Ball
    param[37] = 62.5                    # maximum capacity of Rubisco at 25C-Vcmax	

    param[40] = 0.0015 / RTIMES         # leaf resp co. [kg C-1 d-1 kg-1]
    param[41] = 0.0020 / RTIMES         # stem resp co. [kg C-1 d-1 kg-1]
    param[42] = 0.0020 / RTIMES         # root resp co. [kg C-1 d-1 kg-1]
    param[45] = 0.003 / RTIMES          # fine root resp co. [kg C-1 d-1 kg-1]

    #  the parameter 46 47 were added by G.Mo for the Vcmax - Nitrogen calculations
    param[47] = 3.10 + 1.35   #  leaf Nitrogen content	mean value + 1 SD g/m2 
    param[48] = 20.72 / 62.5   #  slope of Vcmax-N curve  		

  elseif lc == 2  # Set your lc value here for Case 2
    # Conifer Deciduous Forest (DNF) parameters
    param[3] = 0.68                     # clumping_index  
    param[9] = 4.5                      # LAI_max_overstory                     3.3
    param[10] = 2.4                     # LAI_max_understory                    2.4
    param[17] = 1.33                    # z00 roughness length for heat = 0.1 canopy height = 1.5 
    param[19] = 35                      # mass_overstory(kg/m^2)                25
    param[21] = 10                      # mass_understory(kg/m^2)               10
    param[22] = 0.6                     # root_depth(m)                         0.8
    param[23] = 0.035                   # albedo_canopy_vis                     0.09
    param[24] = 0.23                    # albedo_c_nir_observed                 0.19
    param[28] = 0.95                    # decay_rate_of_root_distribution       0.97
    param[29] = 150                     # minimum_stomatal_resistance(s/m)      100
    param[30] = 20                      # canopy_height                         23
    param[31] = 3                       # height of understory                  3

    # the following parameters were added by MEM for the Ball and Berry calculations
    param[34] = 8                       # k Ball
    param[37] = 39.1                    # maximum capacity of Rubisco at 25C-Vcmax	

    param[40] = 0.0015 / RTIMES         # leaf resp co. kg C-1 d-1 kg-1    
    param[41] = 0.0020 / RTIMES         # stem resp co. kg C-1 d-1 kg-1      
    param[42] = 0.0020 / RTIMES         # root resp co. kg C-1 d-1 kg-1   
    param[45] = 0.003 / RTIMES          # fine root resp co. kg C-1 d-1 kg-1	

    # the parameter 46 47 were added by G.Mo for the Vcmax - Nitrogen calculations
    param[47] = 1.81 + 0.64             # leaf Nitrogen content	mean value + 1 SD g/m2 
    param[48] = 22.05 / 39.1            # slope of Vcmax-N curve  		

  elseif lc == 6
    # Broadleaf Deciduous Forest (DBF) parameters
    param[3] = 0.7                      # clumping_index  0.7 
    param[9] = 4.5                      # LAI_max_overstory                     3.3
    param[10] = 2.4                     # LAI_max_understory                    2.4
    param[17] = 1.53                    # z00  roughness length for heat        
    param[19] = 40                      # mass_overstory(kg/m^2)                25
    param[21] = 10                      # mass_understory(kg/m^2)               10
    param[22] = 0.8                     # root_depth(m)                         0.8
    param[23] = 0.04                    # albedo_canopy_vis                     0.09
    param[24] = 0.25                    # albedo_c_nir_observed                 0.19
    param[28] = 0.97                    # decay_rate_of_root_distribution       0.97
    param[29] = 200                     # minimum_stomatal_resistance(s/m)      100
    param[30] = 23                      # canopy_height                         23
    param[31] = 3                       # height of understory                  3

    # the following parameters were added by MEM for the Ball and Berry calculations
    param[34] = 8                       # k Ball			
    param[37] = 57.7                    # maximum capacity of Rubisco at 25C-Vcmax	

    param[40] = 0.015 / RTIMES          # leaf resp co. kg C-1 d-1 kg-1    
    param[41] = 0.0035 / RTIMES         # stem resp co. kg C-1 d-1 kg-1      
    param[42] = 0.0025 / RTIMES         # root resp co. kg C-1 d-1 kg-1   
    param[45] = 0.003 / RTIMES          # fine root resp co.   kg C-1 d-1 kg-1	

    # the parameter 46 47 were added by G.Mo for the Vcmax - Nitrogen calculations
    param[47] = 1.74 + 0.71             # leaf Nitrogen content	mean value + 1 SD g/m2 
    param[48] = 33.79 / 57.7            # slope of Vcmax-N curve  		

  elseif lc == 9
    # Broadleaf Evergreen Forest (EBF) parameters
    param[3] = 0.63                     # clumping_index  0.8 
    param[9] = 4.5                      # LAI_max_overstory                     3.3
    param[10] = 2.4                     # LAI_max_understory                    2.4
    param[17] = 1.53                    # z00  roughness length for heat        
    param[19] = 40                      # mass_overstory(kg/m^2)                25
    param[21] = 10                      # mass_understory(kg/m^2)               10
    param[22] = 0.8                     # root_depth(m)                         0.8
    param[23] = 0.04                    # albedo_canopy_vis                     0.09
    param[24] = 0.25                    # albedo_c_nir_observed                 0.19
    param[28] = 0.97                    # decay_rate_of_root_distribution       0.97
    param[29] = 200                     # minimum_stomatal_resistance(s/m)      100
    param[30] = 23                      # canopy_height                         23
    param[31] = 3                       # height of understory                  3

    # the following parameters were added by MEM for the Ball and Berry calculations
    param[34] = 8                       # k Ball
    param[37] = 29                      # maximum capacity of Rubisco at 25C-Vcmax

    param[40] = 0.015 / RTIMES          # leaf resp co. kg C-1 d-1 kg-1    
    param[41] = 0.0035 / RTIMES         # stem resp co. kg C-1 d-1 kg-1      
    param[42] = 0.0025 / RTIMES         # root resp co. kg C-1 d-1 kg-1   
    param[45] = 0.003 / RTIMES          # fine root resp co.   kg C-1 d-1 kg-1	

    # the parameter 46 47 were added by G.Mo for the Vcmax - Nitrogen calculations
    param[47] = 2.17 + 0.8              # leaf Nitrogen content	mean value + 1 SD g/m2 
    param[48] = 14.02 / 29.0            # slope of Vcmax-N curve  		

  elseif lc == 13
    # Shrub-SH parameters
    param[3] = 0.7                      # clumping_index                        0.8 
    param[9] = 3.3                      # LAI_max_overstory                     3.3
    param[10] = 0.01                    # LAI_max_understory                    2.4
    param[17] = 0.3                     # z00 roughness length for heat         
    param[19] = 25                      # mass_overstory(kg/m^2)                25
    param[21] = 0                       # mass_understory(kg/m^2)               10
    param[22] = 0.5                     # root_depth(m)                         0.8
    param[23] = 0.045                   # albedo_canopy_vis                     0.09
    param[24] = 0.28                    # albedo_c_nir_observed                 0.19
    param[28] = 0.95                    # decay_rate_of_root_distribution       0.97
    param[29] = 100                     # minimum_stomatal_resistance(s/m)      100
    param[30] = 4                       # canopy_height                         23
    param[31] = 0                       # understory_height                     

    # the following parameters were added by MEM for the Ball and Berry calculations
    param[34] = 8                       # k Ball		    	
    param[37] = 61.7 * 0.5 + 54 * 0.5   # maximum capacity of Rubisco at 25C-Vcmax	

    param[40] = 0.001 / RTIMES          # leaf resp co. kg C-1 d-1 kg-1    
    param[41] = 0.002 / RTIMES          # stem resp co. kg C-1 d-1 kg-1      
    param[42] = 0.0015 / RTIMES         # root resp co. kg C-1 d-1 kg-1   
    param[45] = 0.003 / RTIMES          # fine root resp co.   kg C-1 d-1 kg-1

    # the parameter 46 47 were added by G.Mo for the Vcmax - Nitrogen calculations
    param[47] = (2.03 + 1.05 + 1.69 + 0.62) * 0.5   # leaf Nitrogen content	mean value + 1 SD g/m2
    param[48] = (32.09 / 61.7 + 33.14 / 54.0) * 0.5 # slope of Vcmax-N curve

  elseif lc == 40
    # C4 plants parameters
    param[3] = 0.73                     # clumping_index                        0.8 
    param[9] = 4.5                      # LAI_max_overstory                     3.3
    param[10] = 0.01                    # LAI_max_understory                    2.4
    param[17] = 0.04                    # z00 roughness length for heat
    param[19] = 2.5                     # mass_overstory(kg/m^2)                25
    param[21] = 0.1                     # mass_understory(kg/m^2)               10
    param[22] = 0.3                     # root_depth(m)                         0.8
    param[23] = 0.055                   # albedo_canopy_vis                     0.09
    param[24] = 0.3                     # albedo_c_nir_observed                 0.19
    param[28] = 0.95                    # decay_rate_of_root_distribution       0.97
    param[29] = 200                     # minimum_stomatal_resistance(s/m)      100
    param[30] = 4                       # canopy_height                         23
    param[31] = 0.1                     # height of understory                  3

    # the following parameters were added by MEM for the Ball and Berry calculations
    param[34] = 4                       # k Ball			
    param[37] = 30                      # maximum capacity of Rubisco at 25C-Vcmax	

    param[40] = 0.001 / RTIMES          # leaf resp co. kg C-1 d-1 kg-1    
    param[41] = 0.002 / RTIMES          # stem resp co. kg C-1 d-1 kg-1      
    param[42] = 0.0015 / RTIMES         # root resp co. kg C-1 d-1 kg-1   
    param[45] = 0.003 / RTIMES          # fine root resp co.   kg C-1 d-1 kg-1	

    # the parameter 46 47 were added by G.Mo for the Vcmax - Nitrogen calculations
    param[47] = 2.375                   # use crop's leaf Nitrogen content	mean value + 1 SD g/m2 
    param[48] = 0.31                    # use crop's	slope of Vcmax-N curve

  else
    param[3] = 0.8                      # clumping_index                        0.8 
    param[9] = 4.5                      # LAI_max_overstory                     3.3
    param[10] = 2.4                     # LAI_max_understory                    2.4
    param[17] = 0.04                    # z00 roughness length for heat
    param[19] = 2.5                     # mass_overstory(kg/m^2)                25
    param[21] = 0.1                     # mass_understory(kg/m^2)               10
    param[22] = 0.3                     # root_depth(m)                         0.8
    param[23] = 0.055                   # albedo_canopy_vis                     0.09
    param[24] = 0.3                     # albedo_c_nir_observed                 0.19
    param[28] = 0.95                    # decay_rate_of_root_distribution       0.97
    param[29] = 200                     # minimum_stomatal_resistance(s/m)      100
    param[30] = 4                       # canopy_height                         23
    param[31] = 0.1                     # height of understory                  3

    # the following parameters were added by MEM for the Ball and Berry calculations
    param[34] = 8                       # k Ball			
    param[37] = 78.2 * 0.5 + 100.7 * 0.5# maximum capacity of Rubisco at 25C-Vcmax	

    param[40] = 0.001 / RTIMES          # leaf resp co. kg C-1 d-1 kg-1    
    param[41] = 0.002 / RTIMES          # stem resp co. kg C-1 d-1 kg-1      
    param[42] = 0.0015 / RTIMES         # root resp co. kg C-1 d-1 kg-1   
    param[45] = 0.003 / RTIMES          # fine root resp co. kg C-1 d-1 kg-1

    # the parameter 46 47 were added by G.Mo for the Vcmax - Nitrogen calculations
    param[47] = (1.75 + 0.76 + 1.62 + 0.61) * 0.5    # leaf Nitrogen content	mean value + 1 SD g/m2 
    param[48] = (45.29 / 78.2 + 62.75 / 100.7) * 0.5 # slope of Vcmax-N curve
  end

  param[7] = 100                      # light_compensate_point                  30
  param[8] = 1000                     # light_saturation_point                  1000
  param[11] = 0.01                    # LAI_min_overstory                       0.0
  param[12] = 0.01                    # LAI_min_understory                      0.0
  param[13] = 0.87                    # albedo_new_snow                         0.87
  param[14] = 0.94                    # albedo_new_snow_vis                     0.94
  param[15] = 0.8                     # albedo_new_snow_nir                     0.80
  param[16] = 100.0                   # density_now_snow                        100.0
  param[18] = 2700                    # specific_heat_overstory                 2700.0
  param[20] = 2700                    # specific_heat_understory                2700.0
  param[25] = 0.10                    # albedo_saturated_soil                   0.10
  param[26] = 0.35                    # albedo_dry_soil                         0.35
  param[27] = 0.5                     # r_drainage                              1.0
  param[32] = 30                      # the_height_to_measure_wind_speed        30
  param[33] = 0.05                    # the_depth_litter                        0.05
  param[35] = 0.0175                  # intercept_for_H2O_ball_berry
  param[36] = 0.0175 / 1.6            # intercept_for_C2O_ball_berry
  param[38] = 2.39 * param[37] - 14.2 # Jmax at 25 C
  param[39] = 28.0                    # mb / coefficient reflecting the sensitivity of stomata to VPD/moderately N-stressed plants
  param[46] = 2.3                     # Q10 = 2.3 constant for exp. resp.	

  param
end
