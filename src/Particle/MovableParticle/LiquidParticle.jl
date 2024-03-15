#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/09 19:59:58
  @ license: MIT
  @ description:
 =#

mutable struct CommonLiquidParticle{
    RealType <: AbstractFloat,
    ArrayType <: AbstractVector{<:RealType},
    AtomicArrayType <: AbstractVector{Base.Threads.Atomic{RealType}},
} <: LiquidParticle
    x_vec_::ArrayType
    v_vec_::ArrayType
    dv_vec_::AtomicArrayType
    rho_::RealType
    drho_::Base.Threads.Atomic{RealType}
    p_::RealType
    c_::RealType
    mass_::RealType
    gap_::RealType
    sum_kernel_weight_::Base.Threads.Atomic{RealType}
    sum_weighted_scalar_::Base.Threads.Atomic{RealType}
    sum_weighted_vector_::AtomicArrayType
end

function CommonLiquidParticle(dim::IntType where {IntType <: Integer}; RealType::DataType = Float64)::CommonLiquidParticle
    x_vec = zeros(RealType, dim)
    v_vec = zeros(RealType, dim)
    dv_vec = [Base.Threads.Atomic{RealType}(0.0) for _ in 1:dim]
    rho::RealType = 0.0
    drho = Base.Threads.Atomic{RealType}(0.0)
    p::RealType = 0.0
    c::RealType = 0.0
    mass::RealType = 0.0
    gap::RealType = 0.0
    sum_kernel_weight = Base.Threads.Atomic{RealType}(0.0)
    sum_weighted_scalar = Base.Threads.Atomic{RealType}(0.0)
    sum_weighted_vector = [Base.Threads.Atomic{RealType}(0.0) for _ in 1:dim]
    return CommonLiquidParticle(x_vec, v_vec, dv_vec, rho, drho, p, c, mass, gap, sum_kernel_weight, sum_weighted_scalar, sum_weighted_vector)
end

mutable struct ThermalLiquidParticle{
    RealType <: AbstractFloat,
    ArrayType <: AbstractVector{<:RealType},
    AtomicArrayType <: AbstractVector{Base.Threads.Atomic{RealType}},
} <: LiquidParticle
    x_vec_::ArrayType
    v_vec_::ArrayType
    dv_vec_::AtomicArrayType
    rho_::RealType
    drho_::Base.Threads.Atomic{RealType}
    p_::RealType
    c_::RealType
    mass_::RealType
    gap_::RealType
    sum_kernel_weight_::Base.Threads.Atomic{RealType}
    sum_weighted_scalar_::Base.Threads.Atomic{RealType}
    sum_weighted_vector_::AtomicArrayType
    t_::RealType
    dt_::Base.Threads.Atomic{RealType}
    kappa_::RealType
end

function ThermalLiquidParticle(dim::IntType where {IntType <: Integer}; RealType::DataType = Float64)::ThermalLiquidParticle
    x_vec = zeros(RealType, dim)
    v_vec = zeros(RealType, dim)
    dv_vec = [Base.Threads.Atomic{RealType}(0.0) for _ in 1:dim]
    rho::RealType = 0.0
    drho = Base.Threads.Atomic{RealType}(0.0)
    p::RealType = 0.0
    c::RealType = 0.0
    mass::RealType = 0.0
    gap::RealType = 0.0
    sum_kernel_weight = Base.Threads.Atomic{RealType}(0.0)
    sum_weighted_scalar = Base.Threads.Atomic{RealType}(0.0)
    sum_weighted_vector = [Base.Threads.Atomic{RealType}(0.0) for _ in 1:dim]
    t::RealType = 0.0
    dt = Base.Threads.Atomic{RealType}(0.0)
    kappa::RealType = 0.0
    return ThermalLiquidParticle(
        x_vec,
        v_vec,
        dv_vec,
        rho,
        drho,
        p,
        c,
        mass,
        gap,
        sum_kernel_weight,
        sum_weighted_scalar,
        sum_weighted_vector,
        t,
        dt,
        kappa,
    )
end
