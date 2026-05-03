@with_kw mutable struct PhotoConsts{T}
  Γ::T = 0.0
  K::T= 0.0
  Rd_factor::T = 0.0
  Jmax_factor::T = 0.0
  Vcmax_factor::T = 0.0
end

function init_photo_consts(T_leaf_K::T) where {T<:Real}
  tprime25 = T_leaf_K - TK25
  T_factor = tprime25 / (TK25 * rugc * T_leaf_K)

  Kc = kc25 * exp(ekc * T_factor)
  Ko = ko25 * exp(eko * T_factor)
  tau = tau25 * exp(ektau * T_factor)
  Γ = 0.5 * o2 / tau * 1000.0

  K = Kc * (1.0 + o2 / Ko)
  Rd_factor = exp(erd * T_factor)

  Jmax_factor = TBOLTZ(1.0, ejm, toptjm, T_leaf_K)
  Vcmax_factor = TBOLTZ(1.0, evc, toptvc, T_leaf_K)
  Γ, K, Rd_factor, Jmax_factor, Vcmax_factor
end

function PhotoConsts(T_leaf_K::T) where {T<:Real}
  Γ, K, Rd_factor, Jmax_factor, Vcmax_factor = init_photo_consts(T_leaf_K)
  return PhotoConsts{T}(Γ, K, Rd_factor, Jmax_factor, Vcmax_factor)
end

function PhotoConsts!(pc::PhotoConsts{T}, T_leaf_K::T) where {T<:Real}
  Γ, K, Rd_factor, Jmax_factor, Vcmax_factor = init_photo_consts(T_leaf_K)
  @pack! pc = Γ, K, Rd_factor, Jmax_factor, Vcmax_factor
end
