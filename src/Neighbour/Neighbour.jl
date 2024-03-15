#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/09 22:15:17
  @ license: MIT
  @ description:
 =#

abstract type AbstractNeighbour end;

@inline function relativeVelocityVector(
    p_i::ParticleType1 where {ParticleType1 <: MovableParticle},
    p_j::ParticleType2 where {ParticleType2 <: MovableParticle},
)::typeof(p_i.v_vec_)
    return p_i.v_vec_ - p_j.v_vec_
end

@inline function relativeVelocityVector(
    p_i::ParticleType1 where {ParticleType1 <: MovableParticle},
    p_j::ParticleType2 where {ParticleType2 <: FixedParticle},
)::typeof(p_i.v_vec_)
    return p_i.v_vec_
end

@inline function relativeVelocityVector(
    p_i::ParticleType1 where {ParticleType1 <: MovableParticle},
    p_j::ParticleType2 where {ParticleType2 <: StaticVelocityWallParticle},
)::typeof(p_i.v_vec_)
    return p_i.v_vec_ - p_j.v_vec_
end

mutable struct CommonNeighbour{IntType <: Integer, RealType <: AbstractFloat, ArrayType <: AbstractVector{<:RealType}} <: AbstractNeighbour
    i_::IntType
    j_::IntType
    r_::RealType
    r_vec_::ArrayType
    v_vec_::ArrayType
    kernel_value_::RealType
    kernel_gradient_::RealType
    kernel_gradient_vec_::ArrayType
end

function CommonNeighbour(
    neighbour_tuple::Tuple{IntType, IntType, RealType},
    p_i::ParticleType1 where {ParticleType1 <: AbstractParticle},
    p_j::ParticleType2 where {ParticleType2 <: AbstractParticle},
    smooth_kernel::SmoothKernelType where {SmoothKernelType <: SmoothKernel},
)::CommonNeighbour where {IntType <: Integer, RealType <: AbstractFloat}
    i::IntType, j::IntType, r::RealType = neighbour_tuple
    r_vec::typeof(p_i.x_vec_) = p_i.x_vec_ - p_j.x_vec_
    v_vec::typeof(p_i.v_vec_) = relativeVelocityVector(p_i, p_j)
    kernel_value::typeof(r) = kernelValue(r, smooth_kernel)
    kernel_gradient::typeof(r) = kernelGradient(r, smooth_kernel)
    kernel_gradient_vec::typeof(r_vec) = kernel_gradient * normalize(r_vec)
    for i_dim in eachindex(kernel_gradient_vec)
        # avoid nan
        @inbounds if isnan(kernel_gradient_vec[i_dim])
            @inbounds kernel_gradient_vec[i_dim] = 0.0
        end
    end
    return CommonNeighbour(i, j, r, r_vec, v_vec, kernel_value, kernel_gradient, kernel_gradient_vec)
end

function CommonNeighbour!(
    neighbour::CommonNeighbour,
    neighbour_tuple::Tuple{IntType, IntType, RealType},
    p_i::ParticleType1 where {ParticleType1 <: AbstractParticle},
    p_j::ParticleType2 where {ParticleType2 <: AbstractParticle},
    smooth_kernel::SmoothKernelType where {SmoothKernelType <: SmoothKernel},
)::Nothing where {IntType <: Integer, RealType <: AbstractFloat}
    i::IntType, j::IntType, r::RealType = neighbour_tuple
    neighbour.i_ = i
    neighbour.j_ = j
    neighbour.r_ = r
    neighbour.r_vec_ = p_i.x_vec_ - p_j.x_vec_
    neighbour.v_vec_ = relativeVelocityVector(p_i, p_j)
    neighbour.kernel_value_ = kernelValue(r, smooth_kernel)
    neighbour.kernel_gradient_ = kernelGradient(r, smooth_kernel)
    neighbour.kernel_gradient_vec_ = neighbour.kernel_gradient_ * normalize(neighbour.r_vec_)
    for i_dim in eachindex(neighbour.kernel_gradient_vec_)
        # avoid nan
        @inbounds if isnan(neighbour.kernel_gradient_vec_[i_dim])
            @inbounds neighbour.kernel_gradient_vec_[i_dim] = 0.0
        end
    end
    return nothing
end

function findNeighbours!(
    neighbour_tuples::Vector{Tuple{IntType, IntType, RealType}},
    neighbours::Vector{<:AbstractNeighbour},
    ps_i::ParticleArrayType1 where {ParticleArrayType1 <: AbstractVector{<:AbstractParticle}},
    ps_j::ParticleArrayType2 where {ParticleArrayType2 <: AbstractVector{<:AbstractParticle}},
    smooth_kernel::SmoothKernelType where {SmoothKernelType <: SmoothKernel},
)::Nothing where {IntType <: Integer, RealType <: AbstractFloat}
    new_neighbours_number::IntType = length(neighbour_tuples)
    old_neighbours_number::IntType = length(neighbours)
    neighbours = resize!(neighbours, new_neighbours_number)
    if new_neighbours_number <= old_neighbours_number
        @floop @simd for i_neighbour in 1:new_neighbours_number
            @inbounds CommonNeighbour!(
                neighbours[i_neighbour],
                neighbour_tuples[i_neighbour],
                ps_i[neighbour_tuples[i_neighbour][1]],
                ps_j[neighbour_tuples[i_neighbour][2]],
                smooth_kernel,
            )
        end
    else
        @floop @simd for i_neighbour in 1:old_neighbours_number
            @inbounds CommonNeighbour!(
                neighbours[i_neighbour],
                neighbour_tuples[i_neighbour],
                ps_i[neighbour_tuples[i_neighbour][1]],
                ps_j[neighbour_tuples[i_neighbour][2]],
                smooth_kernel,
            )
        end
        @floop @simd for i_neighbour in (old_neighbours_number + 1):new_neighbours_number
            @inbounds neighbours[i_neighbour] = CommonNeighbour(
                neighbour_tuples[i_neighbour],
                ps_i[neighbour_tuples[i_neighbour][1]],
                ps_j[neighbour_tuples[i_neighbour][2]],
                smooth_kernel,
            )
        end
    end
    return nothing
end
