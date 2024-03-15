#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/09 21:13:33
  @ license: MIT
  @ description:
 =#

@inline function updateDensity!(p::ParticleType where {ParticleType <: MovableParticle}, dt::RealType where {RealType <: AbstractFloat})::Nothing
    p.rho_ += p.drho_[] * dt
    p.drho_[] = 0.0
    return nothing
end
