#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/09 21:55:35
  @ license: MIT
  @ description:
 =#

@inline function updatePressure!(
    p::ParticleType where {ParticleType <: FluidParticle},
    wc_lm::WeaklyCompressibleLiquidModelType where {WeaklyCompressibleLiquidModelType <: WeaklyCompressibleLiquidModel},
)::Nothing
    relative_rho::typeof(p.rho_) = p.rho_ / wc_lm.rho_0_
    p.p_ = wc_lm.b_ * (relative_rho^wc_lm.gamma_ - 1) + wc_lm.p_0_
    p.c_ = wc_lm.c_0_ * sqrt(relative_rho^(wc_lm.gamma_ - 1))
    return nothing
end
