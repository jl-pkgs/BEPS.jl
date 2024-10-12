# BEPS.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jl-pkgs.github.io/BEPS.jl/dev)
[![CI](https://github.com/jl-pkgs/BEPS.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/jl-pkgs/BEPS.jl/actions/workflows/CI.yml)
[![Codecov](https://codecov.io/gh/jl-pkgs/BEPS.jl/branch/master/graph/badge.svg)](https://app.codecov.io/gh/jl-pkgs/BEPS.jl/tree/master)

Boreal Ecosystem Productivity Simulator in Julia

> Dongdong Kong
>
> BEPS.jl is alive in Julia now. All functions have been ported to Julia, and the
> performance is about 2.5 times faster than C version.
>
> - Julia: 0.286327 seconds (822.38 k allocations: 22.998 MiB, 0.85% gc time)
> - C    : 0.787059 seconds (629.95 k allocations: 13.915 MiB)

> [!CAUTION]
> `BEPS.clang` only works under Windows.

## Install

- For developers

  ```bash
  git clone https://github.com/jl-pkgs/BEPS.jl
  cd BEPS.jl/deps
  git clone https://github.com/jl-pkgs/BEPS.c
  ```

- For users

  ```bash
  # In Julia
  ] add https://github.com/jl-pkgs/BEPS.jl
  ```

## Usage

```julia
using BEPS
d = deserialize("data/p1_meteo")
lai = readdlm("examples/input/p1_lai.txt")[:]

par = (lon=120.5, lat=30.5, landcover=25, clumping=0.85,
  soil_type=8, Tsoil=2.2,
  soilwater=0.4115, snowdepth=0.0)

@time df_jl, df_ET_jl, Tg = besp_main(d, lai, par; version="julia");
```

> Figure1: The bias of Julia version compared with C, `bias = (Julia - C)/ C`.
![](./docs/images/Figure1_bias_of_julia-version.png)

The bias of `SH` is the largest due to double number accuracy, about 1.48%, which is acceptable.

> Figure2: The variation of soil temperature at different depths.
![](./docs/images/Figure2_variation_of_Tg.png)

See [examples/example_01.qmd](examples/example_01.qmd) for details.

## TODO

- [ ] 通量站上测试模型表现 (top1 task)
- [ ] 热浪期间，土壤温度的表现
- [ ] 土壤类型参数
- [ ] clumping index数据处理

## Bugs

- [ ] snowpack: depth_snow一直增加，一直到无穷

## Researches

- [ ] 研究土壤温度和空气温度之间的关系，为sentinel-2遥感数据反演提供依据
- [ ] 光周期影响测试

## References

1. Hourly BEPS model. <https://github.com/JChen-UToronto/BEPS_hourly_site>

2. Daily BEPS model developed by Jane Liu. <https://github.com/JChen-UToronto/BEPS_D>
