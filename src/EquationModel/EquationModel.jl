#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/09 21:38:19
  @ license: MIT
  @ description:
 =#

abstract type AbstractEquationModel end;

abstract type LiquidModel <: AbstractEquationModel end;
include("./LiquidModel.jl");
