#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/10 22:25:00
  @ license: MIT
  @ description:
 =#

abstract type FixedTimeIntervalDiscretization <: AbstractTimeDiscretization end;

struct DensityReinitializedFixedTimeIntervalForwardEuler{IntType <: Integer, RealType <: AbstractFloat} <: FixedTimeIntervalDiscretization
    dt_::RealType
    total_step_::IntType
    output_interval_::IntType
    check_interval_::IntType
    density_reinitialized_interval_::IntType # recommended to be 30 by SPHysics
end

function isDensityReinitializedStep(
    step::IntType where {IntType <: Integer},
    time_discretization::TimeDiscretizationType where {TimeDiscretizationType <: AbstractTimeDiscretization},
)::Bool
    return step % time_discretization.density_reinitialized_interval_ == 0
end
