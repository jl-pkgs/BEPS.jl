#import "@local/modern-cug-report:0.1.3": *
#show: doc => template(doc, footer: "BEPS.jl", header: "")

// #import "@preview/codly:1.0.0": *// 
// #import "@preview/codedis:0.1.0": code
// #import "@preview/mitex:0.2.5": *
// #import "@preview/showybox:2.0.3": showybox

#show raw: set text(font: ("consolas", "Microsoft Yahei", "SimSun"), size: 10pt)
#set page(numbering: "1")
#set par.line(numbering: "1")
#set par(leading: 0.7em, spacing: 1.24em)
// #set heading(numbering: "1.1")
// #counter(heading).update(0)


#let disp(f) = {
  let dir-root = "../../src/"
  let file = dir-root + f
  let evap = read(file)
  raw(evap, lang: "julia")
  v(1em)
}

= BEPS.jl  // <!-- omit in toc -->

版本：0.1.9，2026-3-25

作者：Dongdong Kong, #link("kongdongdong@cug.edu.cn")


// aerodynamic_conductance.jl
= 1 模型主程序

#disp("BEPS.jl")
#disp("beps_main.jl")
#disp("beps_modern.jl")

= 2 模型核心过程

#disp("BEPS_modules.jl")
#disp("aerodynamic_conductance_V2.jl")
#disp("beps_optimize.jl")
#disp("evaporation_canopy.jl")
#disp("evaporation_soil.jl")

#disp("heat_H_and_LE.jl")
#disp("inter_prg.jl")
#disp("netRadiation.jl")
// #disp("Photosynthesis")
#disp("photosynthesis.jl")
#disp("photosynthesis_helper.jl")
#disp("rainfall_stage.jl")
#disp("snowpack.jl")
#disp("surface_temperature.jl")


= 3 模型数据类型

#disp("DataType/BEPS_State.jl")
#disp("DataType/CanopyLayer.jl")
#disp("DataType/Constant.jl")
#disp("DataType/DataType.jl")
#disp("DataType/Met.jl")
#disp("DataType/OUTPUT.jl")
#disp("DataType/setup.jl")

= 4 模型参数管理

#disp("DataType/Params/BEPS_Param.jl")
#disp("DataType/Params/GlobalData.jl")
#disp("DataType/Params/macro.jl")
#disp("DataType/Params/Param_Init.jl")
#disp("DataType/Params/ParamPhoto.jl")
#disp("DataType/Params/Params.jl")
// data
// deprecated

// clang
// DataType

// SoilPhysics
// SPAC
// standalone
