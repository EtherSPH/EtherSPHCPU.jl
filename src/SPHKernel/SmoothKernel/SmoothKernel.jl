#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/09 17:12:00
  @ license: MIT
  @ description:
 =#

abstract type SmoothKernel <: AbstractSPHKernel end;

include("./CubicSpline.jl");
include("./Gaussian.jl");
include("./WendlandC2.jl");
include("./WendlandC4.jl");

const kSmoothKernelParametersDict = Dict(
    CubicSpline => (2.0, [2.0 / 3.0, 10.0 / 7.0 / pi, 1.0 / pi]),
    Gaussian => (3.0, [1.0 / sqrt(pi), 1.0 / pi, 1.0 / sqrt(pi^3)]),
    WendlandC2 => (2.0, [0.0, 7.0 / 4.0 / pi, 21.0 / 16.0 / pi]),
    WendlandC4 => (2.0, [5.0 / 8.0, 9.0 / 4.0 / pi, 495.0 / 256.0 / pi]),
);

"""
    SmoothKernel(
        influence_radius::RealType where RealType <: AbstractFloat,
        dim::IntType where IntType <: Integer,
        smooth_kernel::UnionAll
    )::smooth_kernel

The `SmoothKernel()` function is used to create a smooth kernel struct.
Where `UnionAll` can be chosen from:

- `CubicSpline`
- `Gaussian`
- `WendlandC2`
- `WendlandC4`
"""
function SmoothKernel(
    influence_radius::RealType where {RealType <: AbstractFloat},
    dim::IntType where {IntType <: Integer},
    smooth_kernel::UnionAll,
)::smooth_kernel
    @assert dim in [1, 2, 3]
    radius_ratio::typeof(influence_radius) = kSmoothKernelParametersDict[smooth_kernel][1]
    h::typeof(influence_radius) = influence_radius / radius_ratio
    sigma::typeof(influence_radius) = kSmoothKernelParametersDict[smooth_kernel][2][dim] / h^dim
    kernel_0::typeof(influence_radius) = sigma
    return smooth_kernel{typeof(dim), typeof(influence_radius)}(h, dim, radius_ratio, influence_radius, sigma, kernel_0)
end
