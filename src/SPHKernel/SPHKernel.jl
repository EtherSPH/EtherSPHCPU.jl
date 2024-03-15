#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/09 17:08:46
  @ license: MIT
  @ description:
 =#

"""
    AbstractSPHKernel

The `AbstractSPHKernel` abstract type is used to define the abstract type of SPH kernel.
It's the parent of all SPH kernel types.
"""
abstract type AbstractSPHKernel end;

include("./SmoothKernel/SmoothKernel.jl");
