#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/10 16:41:49
  @ license: MIT
  @ description:
 =#

@doc raw"""
    reconstructScalar!(
        particles::ParticleArrayType where ParticleArrayType <: AbstractVector{<:AbstractParticle},
        scalar_symbol::Symbol,
        neighbours::NeighbourArrayType where NeighbourArrayType <: AbstractVector{<:AbstractNeighbour},
        smooth_kernel::SmoothKernelType where SmoothKernelType <: SmoothKernel
    )::Nothing

Reconstruct the scalar field of particles using the weighted average of the scalar field of neighbours.

```math
\bar{f}_i = \sum_j \frac{m_j}{\rho_j} f_j W_{ij} / \sum_j \frac{m_j}{\rho_j} W_{ij}
```
"""
function reconstructScalar!(
    particles::ParticleArrayType where {ParticleArrayType <: AbstractVector{<:AbstractParticle}},
    scalar_symbol::Symbol,
    neighbours::NeighbourArrayType where {NeighbourArrayType <: AbstractVector{<:AbstractNeighbour}},
    smooth_kernel::SmoothKernelType where {SmoothKernelType <: SmoothKernel},
)::Nothing
    @floop @simd for neighbour in neighbours
        @inbounds Threads.atomic_add!(
            particles[neighbour.i_].sum_kernel_weight_,
            neighbour.kernel_value_ * particles[neighbour.j_].mass_ / particles[neighbour.j_].rho_,
        )
        @inbounds Threads.atomic_add!(
            particles[neighbour.i_].sum_weighted_scalar_,
            neighbour.kernel_value_ * particles[neighbour.j_].mass_ / particles[neighbour.j_].rho_ * getfield(particles[neighbour.j_], scalar_symbol),
        )
        @inbounds Threads.atomic_add!(
            particles[neighbour.j_].sum_kernel_weight_,
            neighbour.kernel_value_ * particles[neighbour.i_].mass_ / particles[neighbour.i_].rho_,
        )
        @inbounds Threads.atomic_add!(
            particles[neighbour.j_].sum_weighted_scalar_,
            neighbour.kernel_value_ * particles[neighbour.i_].mass_ / particles[neighbour.i_].rho_ * getfield(particles[neighbour.i_], scalar_symbol),
        )
    end
    @floop @simd for particle in particles
        Threads.atomic_add!(particle.sum_kernel_weight_, smooth_kernel.kernel_0_ * particle.mass_ / particle.rho_)
        Threads.atomic_add!(
            particle.sum_weighted_scalar_,
            smooth_kernel.kernel_0_ * particle.mass_ / particle.rho_ * getfield(particle, scalar_symbol),
        )
        setfield!(particle, scalar_symbol, particle.sum_weighted_scalar_[] / particle.sum_kernel_weight_[])
        particle.sum_kernel_weight_[] = 0.0
        particle.sum_weighted_scalar_[] = 0.0
    end
    return nothing
end
