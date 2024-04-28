#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/11 21:16:10
  @ license: MIT
  @ description:
 =#

defaultCheck(p::ParticleType where {ParticleType <: AbstractParticle})::Bool = true;
defaultUserDefinedFunction(; t::RealType where {RealType <: AbstractFloat} = 0.0)::Nothing = nothing;
defaultModify(p::ParticleType where {ParticleType <: AbstractParticle}; t::RealType where {RealType <: AbstractFloat} = 0.0)::Nothing = nothing;

include("./FluidWallSolve.jl");
include("./FluidThermostaticWallSolve.jl");
include("./FluidCompulsiveThermostaticWallSolve.jl");
include("./FluidVelocitySolve.jl");
include("./FluidInletOutletWallSolve.jl");
