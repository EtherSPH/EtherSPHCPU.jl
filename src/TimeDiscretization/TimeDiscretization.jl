#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/10 22:18:42
  @ license: MIT
  @ description:
 =#

abstract type AbstractTimeDiscretization end;

include("./FixedTimeIntervalDiscretization.jl");
abstract type VariableTimeIntervalDiscretization <: AbstractTimeDiscretization end;

function isOutputStep(
    step::IntType where {IntType <: Integer},
    time_discretization::TimeDiscretizationType where {TimeDiscretizationType <: AbstractTimeDiscretization},
)::Bool
    return step % time_discretization.output_interval_ == 0
end

function isCheckStep(
    step::IntType where {IntType <: Integer},
    time_discretization::TimeDiscretizationType where {TimeDiscretizationType <: AbstractTimeDiscretization},
)::Bool
    return step % time_discretization.check_interval_ == 0
end
