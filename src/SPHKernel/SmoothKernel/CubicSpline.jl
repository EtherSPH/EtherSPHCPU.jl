#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/09 17:14:07
  @ license: MIT
  @ description:
 =#

struct CubicSpline{IntType <: Integer, RealType <: AbstractFloat} <: SmoothKernel
    h_::RealType
    dim_::IntType
    radius_ratio_::RealType
    influence_radius_::RealType
    sigma_::RealType
    kernel_0_::RealType
end

@doc raw"""
    kernelValue(
        r::RealType where RealType <: AbstractFloat,
        kernel::CubicSpline
    )::RealType where RealType <: AbstractFloat

The `kernelValue()` function is used to calculate the kernel value at the distance `r` with the `CubicSpline` kernel.

```math
\begin{aligned}
    W(r, h) = \frac{\sigma}{4}
    \begin{cases}
        \begin{array}{ll}
            3q^2(2q-3) + 1, & 0 \leq q < 1, \\
            (2-q)^3, & 1 \leq q < 2, \\
            0, & q \geq 2,
        \end{array}
    \end{cases}
    \quad \text{where} \quad q = \frac{r}{h}.
\end{aligned}
```
"""
function kernelValue(r::RealType, kernel::CubicSpline)::RealType where {RealType <: AbstractFloat}
    q::RealType = r / kernel.h_
    if q < 1.0
        return kernel.sigma_ / 4.0 * (3.0 * q^2 * (q - 2.0) + 4.0)
    elseif q < 2.0
        return kernel.sigma_ / 4.0 * (2.0 - q)^3
    else
        return 0.0
    end
end

@doc raw"""
    kernelGradient(
        r::RealType where RealType <: AbstractFloat,
        kernel::CubicSpline
    )::RealType where RealType <: AbstractFloat

The `kernelGradient()` function is used to calculate the kernel gradient at the distance `r` with the `CubicSpline` kernel.

```math
\begin{aligned}
    \frac{dW}{dr} = \frac{dW}{dq}\frac{dq}{dr} = \frac{1}{h}\frac{dW}{dq}
\end{aligned}
```
"""
function kernelGradient(r::RealType, kernel::CubicSpline)::RealType where {RealType <: AbstractFloat}
    q::RealType = r / kernel.h_
    if q < 1.0
        return kernel.sigma_ / kernel.h_ * 3.0 / 4.0 * q * (3.0 * q - 4.0)
    elseif q < 2.0
        return -kernel.sigma_ / kernel.h_ * 3.0 / 4.0 * (2.0 - q)^2
    else
        return 0.0
    end
end
