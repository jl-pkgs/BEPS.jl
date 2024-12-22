function lai2(Ω::Float64, CosZs::Float64,
  stem_o::Float64, stem_u::Float64,
  lai_o::Float64, lai_u::Float64,
  LAI::Leaf, PAI::Leaf)

  PAI.o_sunlit = CosZs > 0 ? 2 * CosZs * (1 - exp(-0.5 * Ω * (lai_o + stem_o) / CosZs)) : 0
  PAI.o_shaded = (lai_o + stem_o) - PAI.o_sunlit

  PAI.u_sunlit = CosZs > 0 ? 2 * CosZs * (1 - exp(-0.5 * Ω * (lai_o + stem_o + lai_u + stem_u) / CosZs)) - PAI.o_sunlit : 0
  PAI.u_shaded = (lai_u + stem_u) - PAI.u_sunlit

  LAI.o_sunlit = CosZs > 0 ? 2 * CosZs * (1 - exp(-0.5 * Ω * lai_o / CosZs)) : 0
  LAI.o_shaded = max(0, lai_o - LAI.o_sunlit)  # edited by J. Leng

  LAI.u_sunlit = CosZs > 0 ? 2 * CosZs * (1 - exp(-0.5 * Ω * (lai_o + lai_u) / CosZs)) - LAI.o_sunlit : 0
  LAI.u_shaded = max(0, lai_u - LAI.u_sunlit)  # edited by J. Leng
end

function lai2(Ω::Float64, CosZs::Float64,
  stem_o::Float64, stem_u::Float64,
  lai_o::Float64, lai_u::Float64)

  LAI = Leaf()
  PAI = Leaf()
  
  lai2(Ω, CosZs, stem_o, stem_u, lai_o, lai_u, LAI, PAI)
  LAI, PAI
end


function partition_lai(lai, Ω, CosZs)
  sunlit = CosZs > 0 ? 2 * CosZs * (1 - exp(-0.5 * Ω * lai / CosZs)) : 0
  shaded = max(0, lai - sunlit)
  sunlit, shaded
end


"""
    s_coszs(jday::Int, j::Int, lat::Float64, lon::Float64)

# Example
```julia
jday, hour, lat, lon = 20, 12, 20., 120.
s_coszs(jday, hour, lat, lon)
```
"""
function s_coszs(jday::Int, hour::Int, lat::Float64, lon::Float64)
  Delta = 0.006918 - 0.399912 * cos(jday * 2π / 365.0) + 0.070257 * sin(jday * 2π / 365.0) - 
    0.006758 * cos(jday * 4π / 365.0) + 0.000907 * sin(jday * 4π / 365.0)
  # delta is the declination angle of sun.

  hr = hour  + lon / 15.0  # UTC time
  # hr =j*24.0/RTIMES; # local time
  hr > 24 && (hr = hr - 24)
  hr < 0 && (hr = 24 + hr)

  Lat_arc = π * lat / 180.0
  Hsolar1 = (hr - 12.0) * 2.0 * π / 24.0 # local hour angle in arc.

  # sin(h)
  CosZs = cos(Delta) * cos(Lat_arc) * cos(Hsolar1) + sin(Delta) * sin(Lat_arc)
  return CosZs
end
