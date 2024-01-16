# BEPS.jl

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://CUG-hydro.github.io/BEPS.jl/stable) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://CUG-hydro.github.io/BEPS.jl/dev)
[![CI](https://github.com/CUG-hydro/BEPS.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/CUG-hydro/BEPS.jl/actions/workflows/CI.yml)
[![Codecov](https://codecov.io/gh/CUG-hydro/BEPS.jl/branch/master/graph/badge.svg)](https://app.codecov.io/gh/CUG-hydro/BEPS.jl/tree/master)

Boreal Ecosystem Productivity Simulator in Julia

> Dongdong Kong

> [!CAUTION]
> `BEPS.clang` only works under Windows.

## Install
```bash
git clone https://github.com/CUG-hydro/BEPS.jl
cd BEPS.jl/deps
git clone https://github.com/CUG-hydro/BEPS.c
```

## Usage

```julia
using BEPS
d = deserialize("data/p1_meteo")
lai = readdlm("examples/input/p1_lai.txt")[:]

par = (lon=120.5, lat=30.5, landcover=25, clumping=0.85,
  soil_type=8, Tsoil=2.2,
  soilwater=0.4115, snowdepth=0.0)

@time df_jl, df_ET_jl = besp_main(d, lai, par; version="julia");
```

> Figure1: The bias of Julia version compared with C.
![Figure1: bias of julia version compared with c](./images/Figure1_bias_of_julia-version.png)

The bias of `SH` is the largest due to double number accuracy, about 1.48%, which is acceptable.

> The Julia version is about 2X times faster than C.

- Julia:   0.327014 seconds (1.09 M allocations: 47.056 MiB, 2.32% gc time)
- C    :   0.831637 seconds (629.96 k allocations: 13.916 MiB)

See [examples/example_01.qmd](examples/example_01.qmd) for details.

## TODO

- [ ] 土壤类型参数
- [ ] clumping index数据处理
- [ ] 通量站上测试模型表现

## Researches

- [ ] 研究土壤温度和空气温度之间的关系，为sentinel-2遥感数据反演提供依据
- [ ] 光周期影响测试

## References

1. Hourly BEPS model. <https://github.com/JChen-UToronto/BEPS_hourly_site>

2. Daily BEPS model developed by Jane Liu. <https://github.com/JChen-UToronto/BEPS_D>
