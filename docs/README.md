# BEPS.jl 参数管理框架文档

## 📚 核心文档

### [参数管理核心要点](参数管理核心要点.md)
**一页纸搞定** - 所有你需要知道的

## 🚀 快速开始

## 🎯 核心概念（30秒速览）

### 参数定义
```julia
LAI_max_o::FT = 4.5 | (0.1, 7.0)
#       │        │    └─ 边界 (最小值, 最大值)
#       │        └─ 默认值
#       └─ 字段名
```

### 默认值来源
- **硬编码**: `Params.jl` (通用常数)
- **JSON**: `ParamVeg.json` (植被参数)
- **数据表**: `GlobalData.jl` (土壤参数)

### 常用操作
```julia
model = ParamBEPS(VegType=6, SoilType=4; N=5)  # 创建
df = parameters(model)                        # 提取
update!(model, [:veg, :VCmax25], 100.0)      # 更新
```

## 🔗 相关资源

- [AGENTS.md](../AGENTS.md) - AI 智能体指南
- [examples/example_01.qmd](../examples/example_01.qmd) - 基本用法
- [src/DataType/Params/](../src/DataType/Params/) - 源代码
