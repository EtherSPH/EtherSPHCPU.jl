#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/09 19:55:29
  @ license: MIT
  @ description:
 =#

abstract type MovableParticle <: AbstractParticle end;

abstract type FluidParticle <: MovableParticle end;
abstract type LiquidParticle <: FluidParticle end;
include("./LiquidParticle.jl");
# abstract type GasParticle <: FluidParticle end;

# abstract type SolidParticle <: MovableParticle end;
# abstract type RigidParticle <: SolidParticle end;
