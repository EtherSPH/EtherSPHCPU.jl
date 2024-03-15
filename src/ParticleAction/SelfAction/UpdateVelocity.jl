#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/09 21:16:51
  @ license: MIT
  @ description:
 =#

@inline function updateVelocity!(p::ParticleType where {ParticleType <: MovableParticle}, dt::RealType where {RealType <: AbstractFloat})::Nothing
    for i_dim in eachindex(p.v_vec_)
        @inbounds p.v_vec_[i_dim] += p.dv_vec_[i_dim][] * dt
    end
    return nothing
end

@inline function updateVelocity!(
    p::ParticleType where {ParticleType <: MovableParticle},
    dt::RealType,
    body_force_vec::ArrayType where {ArrayType <: AbstractVector{<:RealType}},
)::Nothing where {RealType <: AbstractFloat}
    @inbounds p.v_vec_ .+= body_force_vec .* dt
    updateVelocity!(p, dt)
    return nothing
end
