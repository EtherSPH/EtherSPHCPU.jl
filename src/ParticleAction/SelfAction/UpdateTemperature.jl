#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/09 21:30:24
  @ license: MIT
  @ description:
 =#

@inline function updateTemperature!(p::ParticleType where {ParticleType <: MovableParticle}, dt::RealType where {RealType <: AbstractFloat})::Nothing
    p.t_ += dt * p.dt_[]
    p.dt_[] = 0.0
    return nothing
end
