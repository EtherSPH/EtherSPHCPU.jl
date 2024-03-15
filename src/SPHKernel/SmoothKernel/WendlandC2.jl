#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/09 17:19:18
  @ license: MIT
  @ description:
 =#

struct WendlandC2{IntType <: Integer, RealType <: AbstractFloat} <: SmoothKernel
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
        kernel::WendlandC2
    )::RealType where RealType <: AbstractFloat

The `kernelValue()` function is used to calculate the kernel value at the distance `r` with the `WendlandC2` kernel.

```math
\begin{aligned}
    W(r, h) = \frac{\sigma}{16}
    \begin{cases}
        \begin{array}{ll}
            (2-q)^4(1+2q), & 0 \leq q < 2, \\
            0, & q \geq 2,
        \end{array}
    \end{cases}
    \quad \text{where} \quad q = \frac{r}{h}.
\end{aligned}
```
"""
function kernelValue(r::RealType, kernel::WendlandC2)::RealType where {RealType <: AbstractFloat}
    q::RealType = r / kernel.h_
    if q < 2.0
        return kernel.sigma_ * (2.0 - q)^4 * (1.0 + 2.0 * q) / 16.0
    else
        return 0.0
    end
end

@doc raw"""
    kernelGradient(
        r::RealType where RealType <: AbstractFloat,
        kernel::WendlandC2
    )::RealType where RealType <: AbstractFloat

The `kernelGradient()` function is used to calculate the kernel gradient at the distance `r` with the `WendlandC2` kernel.
See also [kernelGradient](@ref).
"""
function kernelGradient(r::RealType, kernel::WendlandC2)::RealType where {RealType <: AbstractFloat}
    q::RealType = r / kernel.h_
    if q < 2.0
        return -kernel.sigma_ / kernel.h_ * 5.0 / 8.0 * q * (2.0 - q)^3
    else
        return 0.0
    end
end
