const RTIMES = 24.0

"""
```julia
jday = 1
lon = 30
lat = 30
j = 1
s_coszs(jday, j, lat, lon)
```
"""
function s_coszs(jday, j, lat, lon)
  Delta = 0.006918 - 0.399912 * cos(jday * 2.0 * 3.1415926 / 365.0) +
          0.070257 * sin(jday * 2.0 * 3.1415926 / 365.0) -
          0.006758 * cos(jday * 4.0 * 3.1415926 / 365.0) + 0.000907 * sin(jday * 4.0 * 3.1415926 / 365.0)
  # delta is the declination angle of sun.

  hr = j * 24.0 / RTIMES + lon / 15.0  # UTC time
  # //hr =j*24.0/RTIMES; // local time
  hr = clamp(hr, 0, 24)
  # if (hr > 24) hr = 24 - hr; end
  # if (hr < 0) hr = 24 + hr;
  Lat_arc = 3.1415926 * lat / 180.0
  Hsolar1 = (hr - 12.0) * 2.0 * 3.1415926 / 24.0 # local hour angle in arc

  cos(Delta) * cos(Lat_arc) * cos(Hsolar1) + sin(Delta) * sin(Lat_arc)
end

# hPa
"""
saturated pressure

# Arguments

- `t`: temperature, in Kalvin

ES(273.15 + 30)
"""
function ES(t)
  if (t > 0)
    tmp = (54.8781919 - 6790.4985 / t - 5.02808 * log(t))
    es = exp(tmp)
  else
    println("bad es calc")
    es = NaN
  end
  return es
end


