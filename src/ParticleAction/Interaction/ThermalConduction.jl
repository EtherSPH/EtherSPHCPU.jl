#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/10 20:46:01
  @ license: MIT
  @ description:
 =#

@inline function thermalConduction!(
    p_i::ParticleType1 where {ParticleType1 <: AbstractParticle},
    p_j::ParticleType2 where {ParticleType2 <: AbstractParticle},
    neighbour::NeighbourType where {NeighbourType <: AbstractNeighbour},
    smooth_kernel::SmoothKernelType where {SmoothKernelType <: SmoothKernel},
    equation_model::ModelType where {ModelType <: AbstractEquationModel},
)::Nothing
    equivalent_kappa::eltype(p_i.x_vec_) = 2 * p_i.kappa_ * p_j.kappa_ / (p_i.kappa_ + p_j.kappa_)
    dq::eltype(p_i.x_vec_) =
        2 * equivalent_kappa * neighbour.r_ * neighbour.kernel_gradient_ / p_i.rho_ / p_j.rho_ / (neighbour.r_^2 + 0.01 * smooth_kernel.h_^2) *
        (p_i.t_ - p_j.t_)
    Threads.atomic_add!(p_i.dt_, dq * p_j.mass_ / p_i.cp_)
    Threads.atomic_add!(p_j.dt_, -dq * p_i.mass_ / p_j.cp_)
    return nothing
end

@inline function thermalConduction!(
    p_i::ParticleType1 where {ParticleType1 <: FluidParticle},
    p_j::ParticleType2 where {ParticleType2 <: ThermostaticWallParticle},
    neighbour::NeighbourType where {NeighbourType <: AbstractNeighbour},
    smooth_kernel::SmoothKernelType where {SmoothKernelType <: SmoothKernel},
    th_wc_lm::ThermalWeaklyCompressibleLiquidModel,
)::Nothing
    dq::eltype(p_i.x_vec_) =
        1 / th_wc_lm.cp_ * (p_i.kappa_ + p_j.kappa_) * neighbour.r_ * neighbour.kernel_gradient_ / p_i.rho_ / th_wc_lm.rho_0_ /
        (neighbour.r_^2 + 0.01 * smooth_kernel.h_^2) * (p_i.t_ - p_j.t_)
    Threads.atomic_add!(p_i.dt_, dq * p_i.mass_)
    return nothing
end
