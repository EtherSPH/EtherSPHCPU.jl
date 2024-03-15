#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/09 20:26:04
  @ license: MIT
  @ description:
 =#

mutable struct CompulsiveWallParticle{RealType <: AbstractFloat, ArrayType <: AbstractVector{<:RealType}} <: WallParticle
    x_vec_::ArrayType
    normal_vec_::ArrayType
    gap_::RealType
end

function CompulsiveWallParticle(dim::IntType where {IntType <: Integer}; RealType::DataType = Float64)::CompulsiveWallParticle
    x_vec = zeros(RealType, dim)
    normal_vec = zeros(RealType, dim)
    gap::RealType = 0.0
    return CompulsiveWallParticle(x_vec, normal_vec, gap)
end

mutable struct StaticVelocityWallParticle{RealType <: AbstractFloat, ArrayType <: AbstractVector{<:RealType}} <: WallParticle
    x_vec_::ArrayType
    normal_vec_::ArrayType
    gap_::RealType
    v_vec_::ArrayType
end

function StaticVelocityWallParticle(dim::IntType where {IntType <: Integer}; RealType::DataType = Float64)::StaticVelocityWallParticle
    x_vec = zeros(RealType, dim)
    normal_vec = zeros(RealType, dim)
    gap::RealType = 0.0
    v_vec = zeros(RealType, dim)
    return StaticVelocityWallParticle(x_vec, normal_vec, gap, v_vec)
end

mutable struct ThermostaticWallParticle{RealType <: AbstractFloat, ArrayType <: AbstractVector{<:RealType}} <: WallParticle
    x_vec_::ArrayType
    normal_vec_::ArrayType
    gap_::RealType
    t_::RealType
    kappa_::RealType
end

function ThermostaticWallParticle(dim::IntType where {IntType <: Integer}; RealType::DataType = Float64)::ThermostaticWallParticle
    x_vec = zeros(RealType, dim)
    normal_vec = zeros(RealType, dim)
    gap::RealType = 0.0
    t::RealType = 0.0
    kappa::RealType = 0.0
    return ThermostaticWallParticle(x_vec, normal_vec, gap, t, kappa)
end
