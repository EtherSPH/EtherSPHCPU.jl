#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/10 17:40:18
  @ license: MIT
  @ description:
 =#

@inline function pressureForce!(
    p_i::ParticleType1 where {ParticleType1 <: FluidParticle},
    p_j::ParticleType2 where {ParticleType2 <: FluidParticle},
    neighbour::NeighbourType where {NeighbourType <: AbstractNeighbour},
    smooth_kernel::SmoothKernelType where {SmoothKernelType <: SmoothKernel},
)::Nothing
    p_rho_2::typeof(p_i.p_) = p_i.p_ / p_i.rho_^2 + p_j.p_ / p_j.rho_^2
    p_rho_2 += abs(p_rho_2) * 0.01 * neighbour.kernel_value_ / kernelValue((p_i.gap_ + p_j.gap_) / 2, smooth_kernel)
    @simd for i_dim in eachindex(p_i.v_vec_)
        @inbounds Threads.atomic_add!(p_i.dv_vec_[i_dim], -p_j.mass_ * p_rho_2 * neighbour.kernel_gradient_vec_[i_dim])
        @inbounds Threads.atomic_add!(p_j.dv_vec_[i_dim], p_i.mass_ * p_rho_2 * neighbour.kernel_gradient_vec_[i_dim])
    end
    return nothing
end

@inline function pressureForce!(
    p_i::ParticleType where {ParticleType <: FluidParticle},
    neighbour::NeighbourType where {NeighbourType <: AbstractNeighbour},
    smooth_kernel::SmoothKernelType where {SmoothKernelType <: SmoothKernel},
    wc_lm::WeaklyCompressibleModelType where {WeaklyCompressibleModelType <: WeaklyCompressibleLiquidModel},
)::Nothing
    p_rho_2::typeof(p_i.p_) = p_i.p_ / p_i.rho_^2 + wc_lm.p_0_ / wc_lm.rho_0_^2
    p_rho_2 += abs(p_rho_2) * 0.01 * neighbour.kernel_value_ / kernelValue(p_i.gap_, smooth_kernel)
    @simd for i_dim in eachindex(p_i.v_vec_)
        @inbounds Threads.atomic_add!(p_i.dv_vec_[i_dim], -p_i.mass_ * p_rho_2 * neighbour.kernel_gradient_vec_[i_dim])
    end
    return nothing
end
