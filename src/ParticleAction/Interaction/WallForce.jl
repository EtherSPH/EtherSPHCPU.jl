#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/10 18:13:39
  @ license: MIT
  @ description:
 =#

@inline function wallForce!(
    p_i::ParticleType1 where {ParticleType1 <: FluidParticle},
    p_j::ParticleType2 where {ParticleType2 <: WallParticle},
    neighbour::NeighbourType where {NeighbourType <: AbstractNeighbour},
    smooth_kernel::SmoothKernelType where {SmoothKernelType <: SmoothKernel},
)::Nothing
    # Roger & Dalrmple, 2008
    psi::eltype(p_i.x_vec_) = abs(dot(neighbour.r_vec_, p_j.normal_vec_))
    xi::eltype(p_i.x_vec_) = sqrt(max(0, neighbour.r_^2 - psi^2))
    q::eltype(p_i.x_vec_) = psi / p_j.gap_
    if q > 1 || xi > p_j.gap_
        return nothing
    end
    p_xi::eltype(p_i.x_vec_) = abs((1 + cos(pi * xi / p_j.gap_)) / 2)
    verticle_v::eltype(p_i.x_vec_) = dot(neighbour.v_vec_, p_j.normal_vec_)
    beta::eltype(p_i.x_vec_) = verticle_v > 0 ? 0.0 : 1.0
    r_psi::eltype(p_i.x_vec_) = (0.01 * p_i.c_^2 + beta * p_i.c_ * abs(verticle_v)) / smooth_kernel.h_ * abs(1.0 - q) / sqrt(q)
    wall_force::eltype(p_i.x_vec_) = r_psi * p_xi
    @simd for i_dim in eachindex(p_i.v_vec_)
        @inbounds Threads.atomic_add!(p_i.dv_vec_[i_dim], wall_force * p_j.normal_vec_[i_dim])
    end
    return nothing
end
