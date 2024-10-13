"""
    s_coszs(jday::Int, j::Int, lat::Float64, lon::Float64)

# Example
```julia
jday, hour, lat, lon = 20, 12, 20., 120.
s_coszs(jday, hour, lat, lon)
```
"""
function s_coszs(jday::Int, hour::Int, lat::Float64, lon::Float64)
  Delta = 0.006918 - 0.399912 * cos(jday * 2.0 * π / 365.0) + 0.070257 * sin(jday * 2.0 * π / 365.0) - 0.006758 * cos(jday * 4.0 * π / 365.0) + 0.000907 * sin(jday * 4.0 * π / 365.0)
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
