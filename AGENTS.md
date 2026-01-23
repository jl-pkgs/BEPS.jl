# BEPS.jl - Agent Guide

## 0. CRITICAL INSTRUCTION

- **Language:** 请使用中文进行对话

- **Code:** Use English for variable names, comments, and docstrings.

- 代码编写遵循Linux极简主义；同时遵循代码规范、排版美观

This repository contains the Julia implementation of the Boreal Ecosystem Productivity Simulator (BEPS).
It is a land surface model simulating vegetation photosynthesis, energy balance, and soil processes.

## 1. Environment & Setup

- **Language:** Julia 1.10, 1.11
- **Main Branch:** `master`
- **OS Support:** Windows is the primary development environment (BEPS.clang dependency).
- **Performance:** ~2.5x faster than the legacy C version.

## 2. Development Commands

### Setup
```bash
git clone https://github.com/jl-pkgs/BEPS.jl
cd BEPS.jl
# For comparison testing (C version)
cd deps && git clone https://github.com/jl-pkgs/BEPS.c
```

### Testing
```bash
# Activate project
julia --project

# Run all tests
julia --project -e "using Pkg; Pkg.test()"

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

## 3. Code Style & Conventions

**Adhere to "Linux minimalism": Clean, efficient, and beautiful code layout.**

### Naming
- **Greek Letters:** Used extensively for physics (`θ` moisture, `ψ` potential, `κ` conductivity, `ρ` density).
- **Suffixes:** `_p` (previous step), `_c` (current step), `_sat` (saturated), `_vwp` (wilting point).
- **Prefixes:** `f_` (fraction), `r_` (rate), `z_` (depth).
- **Temperature:** Use `T` or `Tsoil`. Do NOT use `temp`.

### Julia Patterns
- **Structs:** Use `@with_kw` for defaults.
- **Fields:** Use `@unpack` to extract and `@pack!` to assign back to structs.
- **Loops:** Use `@inbounds` for performance-critical sections.
- **Math:** Use `@fastmath` where appropriate.
- **Arrays:** Use `.=` for in-place operations to reduce allocations.

### Safety
- **Clamping:** Always clamp physical values to safe ranges.
  - Temperature: `clamp(T, -50, 50)`
  - Moisture: `clamp(θ, θ_vwp, θ_sat)`
  - Ice Ratio: `clamp(ice_ratio, 0.0, 1.0)`

### Documentation
- **Docstrings:** English.
- **Inline Comments:** Chinese is acceptable/encouraged for complex logic explanation.
- **References:** Cite papers (e.g., "Campbell 1974").

## 4. Architecture & Data Structures

### Key Structs
- **Soil (`src/Param/Soil.jl`):** Holds state (`θ`, `Tsoil_c`, `ice_ratio`) and parameters (`Ksat`, `ψ_sat`).
- **TransientCache (`src/DataType/TransientCache.jl`):** 
  - **Optimized:** Only stores variables needing history (k-1) access (e.g., `Cs`, `G`).
  - **Local:** Single-timestep values (ET, simple temp vars) are kept as local variables, NOT cached.

### Integration Loop (`src/inter_prg.jl`)
1. **Soil Water Stress:** `soil_water_factor_v2`
2. **Canopy Processes:** Photosynthesis, radiation, ET.
3. **Thermal Properties:** `UpdateSoilThermalConductivity`, `Update_Cs`.
4. **Surface Temp:** `surface_temperature_jl` (Couples Atmosphere → Soil).
5. **Heat Flux:** `UpdateHeatFlux` (Heat conduction).
6. **Root Uptake:** `Root_Water_Uptake`.
7. **Moisture Update:** `UpdateSoilMoisture`.

## 5. Testing Strategy

- **Framework:** `Test` (standard).
- **Methodology:** Output is validated against the C implementation (`BEPS.c`).
- **Unit Tests:** `test/modules/test-Soil.jl`
- **Integration:** `test/test-beps_main.jl`
- **Key Checks:** Mass balance, Energy balance, Physical bounds.

## 6. Language & Communication

- **User Interaction:** **Must use Chinese** for all conversations.
- **Code:** English naming and docstrings.
