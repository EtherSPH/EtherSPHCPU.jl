#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/10 17:34:06
  @ license: MIT
  @ description:
 =#

@inline function continuity!(
    p_i::ParticleType1 where {ParticleType1 <: FluidParticle},
    p_j::ParticleType2 where {ParticleType2 <: FluidParticle},
    neighbour::NeighbourType where {NeighbourType <: AbstractNeighbour},
)::Nothing
    drho::typeof(p_i.rho_) = dot(neighbour.v_vec_, neighbour.kernel_gradient_vec_)
    Threads.atomic_add!(p_i.drho_, drho * p_j.mass_)
    Threads.atomic_add!(p_j.drho_, drho * p_i.mass_)
    return nothing
end
