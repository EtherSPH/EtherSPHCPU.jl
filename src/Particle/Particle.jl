#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/09 19:49:51
  @ license: MIT
  @ description:
 =#

abstract type AbstractParticle end;

include("./MovableParticle/MovableParticle.jl");
include("./FixedParticle/FixedParticle.jl");
