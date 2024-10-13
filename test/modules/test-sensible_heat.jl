@testset "sensible_heat_jl" begin
  T_leaf = Leaf(20.0)
  T_ground = 15.0
  T_air = 22.0
  rh_air = 80.0

  Gheat = Leaf(100.0)
  Gheat_g = 120.0
  lai = Leaf(2.0)

  r1 = clang.sensible_heat_c(T_leaf, T_ground, T_air, rh_air,
    Gheat, Gheat_g, lai)

  r2 = sensible_heat_jl(T_leaf, T_ground, T_air, rh_air,
    Gheat, Gheat_g, lai)
  all(r1 .≈ r2)
end
