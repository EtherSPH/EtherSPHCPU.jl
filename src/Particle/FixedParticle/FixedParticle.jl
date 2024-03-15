#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/09 19:56:16
  @ license: MIT
  @ description:
 =#

abstract type FixedParticle <: AbstractParticle end;

abstract type WallParticle <: FixedParticle end;
include("./WallParticle.jl");
