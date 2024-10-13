@with_kw mutable struct State{FT}
  "Surface Temperature: [Ts0, Tsn, Tsm0, Tsn1, Tsn2]"
  Ts::Vector{FT} = zeros(FT, 5)         # 4:8
  Tsoil_prev::Vector{FT} = zeros(FT, 6) # 10:15
  θ_prev::Vector{FT} = zeros(FT,6)      # 22:27
  ice_ratio::Vector{FT} = zeros(FT,6)   # 28:33

  Qhc_o::FT=0.0
  m_water::Layer2 = Layer2{FT}()     # [15, 18] + 1
  m_snow::Layer3 = Layer3{FT}() # [16, 19, 20] + 1
end


function update_state!(state::State, var_n::Vector{Float64})
  Tsoil = state.Tsoil
  var_n[3+1] = Tsoil[1]       # Ts0, 3
  var_n[4+1] = Tsoil[2]       # Tsn, 4 
  var_n[5+1] = Tsoil[3]       # Tsm0, 5
  var_n[6+1] = Tsoil[4]       # Tsn1, 6
  var_n[7+1] = Tsoil[5]       # Tsn1, 7

  var_n[11+1] = state.Qhc_o   # Qhc_o, 11, sensible heat flux

  var_n[15+1] = state.m_water.o
  var_n[18+1] = state.m_water.u

  var_n[16+1] = state.m_snow.o
  var_n[19+1] = state.m_snow.u
  var_n[20+1] = state.m_snow.g
end

export State
