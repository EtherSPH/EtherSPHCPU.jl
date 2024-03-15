#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/10 17:52:58
  @ license: MIT
  @ description:
 =#

@inline function viscosityForce!(
    p_i::ParticleType1 where {ParticleType1 <: FluidParticle},
    p_j::ParticleType2 where {ParticleType2 <: FluidParticle},
    neighbour::NeighbourType where {NeighbourType <: AbstractNeighbour},
    smooth_kernel::SmoothKernelType where {SmoothKernelType <: SmoothKernel},
    wc_lm::WeaklyCompressibleModelType where {WeaklyCompressibleModelType <: WeaklyCompressibleLiquidModel},
)::Nothing
    sum_rho::typeof(p_i.rho_) = p_i.rho_ + p_j.rho_
    mu_i::typeof(p_i.rho_) = wc_lm.nu_0_ * p_i.rho_
    mu_j::typeof(p_i.rho_) = wc_lm.nu_0_ * p_j.rho_
    viscosity_force::typeof(p_i.rho_) = # r ⋅ ∇W = |r|W'
        4 * (mu_i + mu_j) * neighbour.r_ * neighbour.kernel_gradient_ / sum_rho^2 / (neighbour.r_^2 + 0.01 * smooth_kernel.h_^2)
    @simd for i_dim in eachindex(p_i.v_vec_)
        @inbounds Threads.atomic_add!(p_i.dv_vec_[i_dim], p_j.mass_ * viscosity_force * neighbour.v_vec_[i_dim])
        @inbounds Threads.atomic_add!(p_j.dv_vec_[i_dim], -p_i.mass_ * viscosity_force * neighbour.v_vec_[i_dim])
    end
    return nothing
end

@inline function viscosityForce!(
    p_i::ParticleType1 where {ParticleType1 <: FluidParticle},
    p_j::ParticleType2 where {ParticleType2 <: WallParticle},
    neighbour::NeighbourType where {NeighbourType <: AbstractNeighbour},
    smooth_kernel::SmoothKernelType where {SmoothKernelType <: SmoothKernel},
    wc_lm::WeaklyCompressibleModelType where {WeaklyCompressibleModelType <: WeaklyCompressibleLiquidModel},
)::Nothing
    sum_rho::typeof(p_i.rho_) = 2 * p_i.rho_
    mu_i::typeof(p_i.rho_) = wc_lm.nu_0_ * p_i.rho_
    viscosity_force::typeof(p_i.rho_) = # r ⋅ ∇W = |r|W'
        8 * mu_i * neighbour.r_ * neighbour.kernel_gradient_ / sum_rho^2 / (neighbour.r_^2 + 0.01 * smooth_kernel.h_^2)
    verticle_v::typeof(p_i.rho_) = dot(neighbour.v_vec_, p_j.normal_vec_)
    @simd for i_dim in eachindex(p_i.v_vec_)
        @inbounds Threads.atomic_add!(
            p_i.dv_vec_[i_dim],
            p_i.mass_ * viscosity_force * (neighbour.v_vec_[i_dim] - verticle_v * p_j.normal_vec_[i_dim]),
        )
    end
    return nothing
end
