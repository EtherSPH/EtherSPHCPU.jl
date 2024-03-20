#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/11 21:16:10
  @ license: MIT
  @ description:
 =#

defaultCheck(p::ParticleType where {ParticleType <: AbstractParticle})::Bool = true;

include("./FluidWallSolve.jl");
include("./FluidCompulsiveThermostaticWallSolve.jl");
include("./FluidVelocitySolve.jl");
