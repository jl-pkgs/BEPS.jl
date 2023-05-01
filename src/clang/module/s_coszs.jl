"""
s_coszs(jday::Int, j::Int, lat::T, lon::T) where {T<:Real}

# Example
```jldoc
x = zeros(4)
s_coszs(Ref(x, 1), 10, 1, 30, 120)
x
```
"""
function s_coszs(CosZs,
  jday::Int, j::Int, lat::T, lon::T) where {T<:Real}
  ## This version works
  # CosZs = [0.0]
  # CosZs = init_dbl()
  ccall((:s_coszs, libbeps), Cvoid, (Cshort, Cshort, Cfloat, Cfloat, Ptr{Cdouble}),
    jday, j, lat, lon, CosZs)
  # return Value(CosZs)
end

function s_coszs(jday::Int, j::Int, lat::T, lon::T) where {T<:Real}
  ## This version works
  # CosZs = [0.0]
  CosZs = init_dbl()
  ccall((:s_coszs, libbeps), Cvoid, (Cshort, Cshort, Cfloat, Cfloat, Ptr{Cdouble}),
    jday, j, lat, lon, CosZs)
  # return Value(CosZs)
  CosZs[]
end
