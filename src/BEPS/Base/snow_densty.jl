# LoTmpDnsTruncatedAnderson1976
"""
    snow_density(Ta::Float64, U10::Float64=NaN; tfrz=0.0, method="LoTmpDnsSlater2017")

Calculate the density of new snow

"Snow fraction depends on the density of snow in CLM4. A 10 cm snowpack has f_snow=1 when density is low
(50–100 kg m–3), such as may be found in fresh snow, and a smaller snow fraction (f_snow = 0.76)
when density is high (400 kg m–3)" -- Bonan 2019, P148

Snow compacts over time, increasing to a density of 100–500 kg m–3.

# Reference

- van Kampenhout et al. 2017

- CLM5, <https://github.com/ESCOMP/CTSM/blob/07051e3758addf2f9753d520823be9ebcbfec0aa/src/biogeophys/SnowHydrologyMod.F90#L3729-L3747>

# Examples
```julia
Ta = -100.:50
ρ_snow = snow_density.(Ta, 2.0)
ρ_snow_chen = @. 67.9 + 51.3 * exp(Ta / 2.6)

# using Plots
# plot(Ta, ρ_snow)
# plot!(Ta, ρ_snow_chen)
```
"""
function cal_snow_density(Ta::Float64, U10::Float64=NaN; tfrz=0.0, method="LoTmpDnsSlater2017")
  if Ta > tfrz + 2.0
    ρ_snow = 50.0 + 1.7 * (17.0)^1.5
  elseif Ta > tfrz - 15.0
    ρ_snow = 50.0 + 1.7 * (Ta - tfrz + 15.0)^1.5
  else
    if method == "LoTmpDnsTruncatedAnderson1976"
      ρ_snow = 50.0
    elseif method == "LoTmpDnsSlater2017"
      # ρ_snow = -3.833 * (Ta - tfrz) - 0.0333 * (Ta - tfrz)^2
      t_for_bifall_degC = Ta > tfrz - 57.55 ? (Ta - tfrz) : -57.55
      ρ_snow = -(50.0 / 15.0 + 0.0333 * 15.0) * t_for_bifall_degC - 0.0333 * t_for_bifall_degC^2
    end
  end

  if U10 > 0.1
    ρ_snow = ρ_snow + (266.861 * ((1.0 + tanh(U10 / 2)) / 2.0)^8.8)
  end
  ρ_snow
end
