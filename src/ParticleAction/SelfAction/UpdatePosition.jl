#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/09 21:26:13
  @ license: MIT
  @ description:
 =#

@inline function updatePosition!(p::ParticleType where {ParticleType <: MovableParticle}, dt::RealType where {RealType <: AbstractFloat})::Nothing
    for i_dim in eachindex(p.x_vec_)
        @inbounds p.x_vec_[i_dim] += p.v_vec_[i_dim] * dt - p.dv_vec_[i_dim][] * dt^2 / 2
        @inbounds p.dv_vec_[i_dim][] = 0.0
    end
    return nothing
end
