# BEPS.jl

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://CUG-hydro.github.io/BEPS.jl/stable) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://CUG-hydro.github.io/BEPS.jl/dev)
[![CI](https://github.com/CUG-hydro/BEPS.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/CUG-hydro/BEPS.jl/actions/workflows/CI.yml)
[![Codecov](https://codecov.io/gh/CUG-hydro/BEPS.jl/branch/master/graph/badge.svg)](https://app.codecov.io/gh/CUG-hydro/BEPS.jl/tree/master)

Boreal Ecosystem Productivity Simulator in Julia

> Dongdong Kong

> [!CAUTION]
> `BEPS.clang` only works under Windows.

## Usage
```bash
git clone https://github.com/CUG-hydro/BEPS.jl
cd BEPS.jl/deps
git clone https://github.com/CUG-hydro/BEPS.c
```

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
