"""
    prepare_and_compute_surface_temperature!(soil, var, k, T_air, RH,
        z_snow, z_water, Gheat_g, ρ_snow, Tc_u, radiation_g,
        Evap_soil, Evap_SW, Evap_SS, f_snow_g, dz, κ)

辅助函数：准备数据并计算表面温度

# 功能说明
这个函数封装了以下操作：
1. 从 soil 提取并准备热力学参数（Cs, κ, dz, T, G）
2. 保存到 var 的相应位置
3. 调用 surface_temperature_jl 计算表面温度
4. 将结果保存回 var 和 soil

# 参数
- `soil`: 土壤状态结构体
- `var`: 中间临时变量结构体
- `k`: 当前时间步索引
- 其他参数同 surface_temperature_jl

# 副作用
- 修改 var.Cs, var.Tc_u, var.T_soil, var.G
- 修改 var.G[1,k], var.T_ground[k], var.T_soil[1,k], var.T_surf_mix[k],
  var.T_surf_snow[k], var.T_snow_L1[k], var.T_snow_L2[k]
- 修改 soil.Tsoil_c[1]
- 修改 dz[2], κ[2]（临时数组）
"""
function prepare_and_compute_surface_temperature!(
  soil::AbstractSoil, var::InterTempVars, k::Int,
  T_air::FT, RH::FT, z_snow::FT, z_water::FT,
  Gheat_g::FT, ρ_snow::FT, Tc_u::FT,
  radiation_g::FT, Evap_soil::FT, Evap_SW::FT, Evap_SS::FT,
  f_snow_g::FT, dz::Vector{FT}, κ::Vector{FT}) where {FT<:AbstractFloat}

  # 1. 准备土壤热力学参数
  var.Cs[1:2, k] .= soil.Cs[1]      # 体积热容 [J m-3 K-1]
  var.Tc_u[k] = Tc_u                # 下层冠层温度 [°C]
  κ[2] = soil.κ[1]                   # 热导率 [W m-1 K-1]
  dz[2] = soil.dz[1]                 # 层厚度 [m]

  # 2. 准备上一时刻的温度状态
  var.T_soil[1, k-1] = soil.Tsoil_p[1]  # 第1层土壤温度 [°C]
  var.T_soil[2, k-1] = soil.Tsoil_p[2]  # 第2层土壤温度 [°C]
  var.G[2, k] = soil.G[1]                # 第1层土壤热通量 [W m-2]

  # 3. 调用 surface_temperature_jl 计算
  var.G[1, k], var.T_ground[k], var.T_soil[1, k], var.T_surf_mix[k],
  var.T_surf_snow[k], var.T_snow_L1[k], var.T_snow_L2[k] =
    surface_temperature_jl(T_air, RH, z_snow, z_water,
      var.Cs[2, k], var.Cs[1, k], Gheat_g, dz[2], ρ_snow, var.Tc_u[k],
      radiation_g, Evap_soil, Evap_SW, Evap_SS,
      κ[2], f_snow_g, var.G[2, k],
      var.T_ground[k-1],
      var.T_soil[2, k-1], var.T_soil[1, k-1], var.T_surf_mix[k-1],
      var.T_surf_snow[k-1], var.T_snow_L1[k-1], var.T_snow_L2[k-1])

  # 4. 更新 soil 的当前温度状态
  soil.Tsoil_c[1] = var.T_soil[1, k]

  return nothing
end
