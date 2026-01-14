# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

BEPS (Boreal Ecosystem Productivity Simulator) - Julia implementation of a land surface model simulating vegetation photosynthesis, energy balance, and soil processes.

**Performance:** ~2.5x faster than C version
**Main Branch:** `master`
**Julia Version:** 1.8, 1.9, 1.10

## Development Commands

### Setup
```bash
# Clone repository
git clone https://github.com/jl-pkgs/BEPS.jl
cd BEPS.jl

# For developers - also clone C version for comparison
cd deps
git clone https://github.com/jl-pkgs/BEPS.c
```

### Testing
```bash
# Activate project environment
julia --project

# Run all tests
julia --project -e "using Pkg; Pkg.test()"

# Compile check only
julia --project -e "using BEPS"

# Run specific example
julia --project examples/example_01.qmd
```

### Running the Model
```julia
using BEPS
d = deserialize("data/p1_meteo")
lai = readdlm("examples/input/p1_lai.txt")[:]

par = (lon=120.5, lat=30.5, landcover=25, clumping=0.85,
  soil_type=8, Tsoil=2.2, soilwater=0.4115, snowdepth=0.0)

@time df_jl, df_ET_jl, Tg = besp_main(d, lai, par; version="julia")
```

## Soil Water & Heat Module Summary

### Core Physics Components

#### 1. Soil Moisture (`src/Soil/UpdateSoilMoisture.jl`)
- **Model:** Campbell (1974) 5-layer hydraulic model
- **Process:** Darcy flow with gravity + matric potential gradients
- **Key equations:**
  - Water potential: `ψ = ψ_sat * (θ / θ_sat)^(-b)`
  - Hydraulic conductivity: `K = K_sat * (θ / θ_sat)^(2b + 3)`
- **Features:**
  - Adaptive time-stepping (1-360s based on flow velocity)
  - Green-Ampt style infiltration with saturation limit
  - Root water uptake distributed by layer
  - Surface runoff from excess precipitation
  - Ice content affects hydraulic properties (`f_water` factor)

#### 2. Soil Heat (`src/Soil/UpdateHeatFlux.jl`)
- **Model:** 5-layer heat conduction with freeze-thaw
- **Key functions:**
  - `UpdateHeatFlux()` - Heat conduction between layers
  - `UpdateSoilThermalConductivity()` - Johansen method (dry + ice + water components)
  - `Update_Cs()` - Volumetric heat capacity (soil + water + ice + organic matter)
  - `Update_ice_ratio()` - Phase change with latent heat (Lf0 = 334 kJ/kg)
- **Physics:**
  - Heat flux: `G[i] = 2(T[i-1] - T[i]) / (dz[i-1]/κ[i-1] + dz[i]/κ[i])`
  - Temperature update: `T_new = T_old + (G_in - G_out) / (Cs * dz) * Δt`
  - Phase change holds T at 0°C during freezing/melting

#### 3. Surface Temperature (`src/surface_temperature.jl`)
- **Model:** Energy balance coupling atmosphere to soil
- **Regimes:** 3 cases based on snow depth (≤2cm, 2-5cm, >5cm)
- **Energy equation:** `Rn - LE - H - G = 0`
- **Handles:** Snow layers (up to 3), frozen soil surface, ponded water

#### 4. Soil Parameters (`src/Soil/Init_Soil_Parameters.jl`)
- **Database:** 11 soil texture classes (sand → clay)
- **Parameters:** Ksat, porosity, field capacity, wilting point, ψ_sat, b, κ_dry
- **Layers:** Default depths [0.05, 0.10, 0.20, 0.40, 1.25] m
- **Vegetation:** Adjusts ψ_min and alpha for DBF/EBF (VegType 6/9)

### Soil Struct (`src/Param/Soil.jl`)

**State Variables:**
- `θ[1:10]` - Soil moisture [m³/m³]
- `Tsoil_c[1:10]` - Current temperature [°C]
- `ice_ratio[1:10]` - Frozen fraction [0-1]
- `z_water` - Ponded water depth [m]
- `z_snow` - Snow depth [m]
- `r_rain_g` - Rainfall reaching ground [m/s]

**Parameters:**
- `dz[1:10]` - Layer thickness [m]
- `θ_sat, θ_vwp, Ksat, ψ_sat, b` - Hydraulic properties
- `κ_dry, ρ_soil, V_SOM` - Thermal properties
- `r_drainage` - Surface drainage rate
- `r_root_decay` - Root distribution decay

**Derived Variables:**
- `ψ, km, κ, Cs, G` - Matric potential, conductivity, thermal properties, heat flux
- `Ett[1:10]` - ET per layer [m/s]
- `f_soilwater` - Overall water stress factor

### Integration Pattern (`src/inter_prg.jl` lines 158-328)

**Execution Order per Timestep:**
```julia
# 1. Soil water stress factor
soil_water_factor_v2(soil)

# 2. [Canopy processes: radiation, photosynthesis, transpiration, evaporation]
#    → Outputs: Trans_o, Trans_u, Evap_soil, radiation_g, Gheat_g

# 3. Soil thermal properties
UpdateSoilThermalConductivity(soil)
Update_Cs(soil)

# 4. Surface temperature (couples atm → soil)
G[1], Ts0, Tm[1], ... = surface_temperature_jl(
    Ta, RH, z_snow, z_water, Cs, Gheat_g, dz, ρ_snow, Tc_u,
    radiation_g, Evap_soil, Evap_SW, Evap_SS, κ, f_snow, G[2], ...)

soil.Tsoil_c[1] = Tm[1]
soil.G[1] = G[1]

# 5. Soil heat flux (updates layers 2-5)
UpdateHeatFlux(soil, Ta_annual, kstep)

# 6. Root water uptake distribution
Root_Water_Uptake(soil, Trans_o, Trans_u, Evap_soil)
# → Sets soil.Ett[1:5]

# 7. Soil moisture update (all layers)
UpdateSoilMoisture(soil, kstep)
```

**Key Constants:**
- `kstep = 360.0` seconds (sub-hourly timestep for soil)
- `step = 3600.0` seconds (hourly output)
- `layer = 5` (soil layers)
- `ρ_w = 1025.0` kg/m³
- `ρₐ = 1.292` kg/m³

### Dependencies to Vegetation Modules

**External Inputs:**
1. **Transpiration:** `Trans_o, Trans_u` from photosynthesis module
2. **Soil Evaporation:** `Evap_soil` from `evaporation_soil_jl()`
3. **Net Radiation:** `radiation_g` from `netRadiation_jl()`
4. **Aerodynamic Conductance:** `Gheat_g = 1 / ra_g`
5. **Snow/Water State:** `z_snow, z_water` from snowpack/rainfall modules
6. **Canopy Temperature:** `Tc.u` for surface temperature regime

**Can be Isolated by:**
- Setting `Trans_o = Trans_u = 0` (no plant uptake)
- Using simplified bare soil energy balance for `radiation_g, Gheat_g`
- Simple snow accumulation/melt instead of full snowpack model
- Setting `Tc.u = Ta` (no understory canopy)

## Code Style Conventions

**Naming:**
- Greek letters: `θ` (moisture), `ψ` (potential), `κ` (conductivity), `ρ` (density)
- Suffixes: `_p` (previous), `_c` (current), `_sat` (saturated), `_vwp` (wilting point)
- Prefixes: `f_` (fraction/factor), `r_` (rate), `z_` (depth)
- Temperature: `T` or `Tsoil`, not `temp`

**Julia Patterns:**
- `@unpack` for struct field extraction
- `@pack!` for struct field assignment (replaces manual field-by-field assignment)
- `@inbounds` for performance-critical loops
- `@fastmath` for mathematical functions
- `@with_kw` for struct definitions with defaults
- `.=` for in-place array operations

**Safety/Clamping:**
- Temperature: `clamp(T, -50, 50)` or `clamp(T, Ta-25, Ta+25)`
- Soil moisture: `clamp(θ, θ_vwp, θ_sat)`
- Heat flux: `clamp(G, -200, 200)` or `clamp(G, -100, 100)`
- Ice ratio: `clamp(ice_ratio, 0.0, 1.0)`

**Documentation:**
- Docstrings in English
- Inline comments in Chinese
- References to papers (e.g., "Campbell 1974", "Chen 2007, Eq. 18")

## Data Structure Optimizations (Recent Work)

### TransientCache Simplification
The `TransientCache` struct (src/DataType/TransientCache.jl) has been heavily optimized to remove redundant storage:

**Current structure (9 fields only):**
```julia
@with_kw mutable struct TransientCache
  # Temperature state (need k-1 access for thermal calculations)
  Tc_u::Vector{FT}          # Understory canopy temperature
  T_ground::Vector{FT}      # Ground surface temperature
  T_surf_mix::Vector{FT}    # Mixed surface temperature
  T_surf_snow::Vector{FT}   # Snow surface temperature
  T_snow_L1::Vector{FT}     # Snow layer 1 temperature
  T_snow_L2::Vector{FT}     # Snow layer 2 temperature

  # Soil thermal (multi-layer × timestep)
  Cs::Matrix{FT}            # Soil volumetric heat capacity
  G::Matrix{FT}             # Soil layer heat flux

  # Leaf cache (energy balance iteration state)
  leaf_cache::LeafCache
end
```

**Optimization principles:**
1. **ET variables** - Now local variables in `inter_prg.jl`, not cached
2. **Soil temperature** - Use `soil.Tsoil_p` directly, no `cache.T_soil` matrix
3. **Water/snow masses** - Removed 16 unused vector fields (Wcl_*, Wcs_*, Xcl_*, etc.)
4. **Total savings:** ~4.8 KB per instance, 70% reduction in fields

**Key pattern:** Only cache variables that need historical access (k-1) for thermal calculations. All single-timestep values should be local variables.

### State vs Cache vs Parameters
- **State** (`src/Param/Soil.jl`) - Persistent across hours, saved between runs
- **Cache** (`src/DataType/TransientCache.jl`) - Temporary sub-hourly loop storage
- **Parameters** - Model configuration, constant during simulation

**State struct:**
```julia
@with_kw mutable struct State{FT}
  Ts::Vector{FT}         # Surface temperatures [T_ground, T_surf_snow, ...]
  Qhc_o::FT              # Sensible heat flux (needs previous value)
  m_water::Layer2        # Canopy water mass
  m_snow::Layer3         # Canopy/ground snow mass
  ρ_snow::FT             # Snow density
end
```

## File Structure

```
src/
├── beps_main.jl                # Top-level interface (besp_main function)
├── inter_prg.jl                # Main integration loop (hourly timestep)
├── surface_temperature.jl      # Surface energy balance
├── evaporation_soil.jl         # Bare soil evaporation
├── Soil/
│   ├── UpdateSoilMoisture.jl   # Water balance (core physics)
│   ├── UpdateHeatFlux.jl       # Heat conduction (core physics)
│   ├── Init_Soil_Parameters.jl # Soil texture database
│   └── soil_water_factor_v2.jl # Water stress calculation
├── Param/
│   └── Soil.jl                 # Soil and State struct definitions
├── DataType/
│   ├── OUTPUT.jl               # Results and OutputET structs
│   ├── TransientCache.jl       # Sub-hourly cache (optimized)
│   └── Constant.jl             # Physical constants
└── standalone/
    └── UpdateSoilMoisture.jl   # Simplified version (no f_water, ice updates)

test/
├── test-beps_main.jl           # Integration tests
└── modules/
    └── test-Soil.jl            # Soil module unit tests
```

## Known Issues and Bug Fixes

### Critical Fixes (2024-10-13 ~ 2024-10-14)

**Snowpack bugs** - Snow accumulation issues causing unrealistic depths:
- `snowpack_stage1`: Uninitialized `snowrate_o` causing incorrect conditions
- `snowpack_stage3`: `max` should be `min` for `mass_water_frozen`
- **Fix**: Corrected melt/freeze conditions and added 10m depth limit
  ```julia
  # Correct conditions:
  con_melt = Tsnow > 0 && ms_sup > 0
  con_frozen = Tsnow <= 0 && z_water > 0
  ```
- `ρ_snow`: Initial value set to 250 kg/m³ (now in `state.ρ_snow`)

**Photosynthesis bug** - Incorrect latent heat constant:
- `LAMBDA` in photosynthesis: `lambda_ice` was 333 J/kg, should be 333000 J/kg

### Outstanding Issues (2025-10-25)
- [ ] Leaf temperature parameter passing in photosynthesis module needs correction

## Testing

**Unit tests:** `test/modules/test-Soil.jl`
**Integration:** `test/test-beps_main.jl`

**Key validation checks:**
- Mass balance: `Σ(Δθ * dz) = infiltration - runoff - ET - drainage`
- Energy balance: `Rn - LE - H - G ≈ 0`
- Physical bounds: `θ ∈ [θ_vwp, θ_sat]`, `T ∈ [-50, 50]`, `ice ∈ [0, 1]`
- C/Julia comparison: All soil struct fields match between versions

**Test structure:**
- Tests compare Julia vs C implementation results
- Format: "C and Julia: field_name, field_name" for each validated field
- 39 fields validated in UpdateHeatFlux and Init_Soil_var tests

## Quick Reference

**Initialize soil:**
```julia
soil = Soil()
Init_Soil_Parameters(VegType, SoilType, r_root_decay, soil)
Init_Soil_Status(soil, Tsoil0, Tair, θ0, z_snow0)
```

**Update soil (one timestep):**
```julia
# 1. Thermal properties
UpdateSoilThermalConductivity(soil)
Update_Cs(soil)

# 2. Surface temperature
G, Ts, ... = surface_temperature_jl(...)
soil.Tsoil_c[1] = Ts
soil.G[1] = G

# 3. Heat flux
UpdateHeatFlux(soil, Ta_annual, kstep)

# 4. Water uptake
soil.Ett[1] = Evap_soil / ρ_w
soil.Ett[2:5] .= 0.0

# 5. Moisture
UpdateSoilMoisture(soil, kstep)
```

**Access results:**
```julia
θ = soil.θ[1:5]                # Moisture by layer
T = soil.Tsoil_c[1:5]          # Temperature by layer
ice = soil.ice_ratio[1:5]      # Ice fraction by layer
```

## Testing

**Unit tests:** `test/modules/test-Soil.jl`
**Integration:** `test/test-beps_modern.jl`

**Key checks:**
- Mass balance: `Σ(Δθ * dz) = infiltration - runoff - ET - drainage`
- Energy balance: `Rn - LE - H - G ≈ 0`
- Physical bounds: `θ ∈ [θ_vwp, θ_sat]`, `T ∈ [-50, 50]`, `ice ∈ [0, 1]`

## References

**Physics:**
- Campbell (1974) - Soil water retention and hydraulic conductivity
- Chen (2007) Ecological Modelling - Volumetric heat capacity (Eq. 18)
- Bonan (2019) Table 8.2 - Campbell parameters
- He et al. (2017) JGR-B - Soil water stress factor (Eq. 4-5)

**Model:**
- BEPS V2023 - Current implementation version
- CLM3.5 - Reference for zero-flow bottom BC
