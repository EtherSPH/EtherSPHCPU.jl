#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/09 17:20:04
  @ license: MIT
  @ description:
 =#

struct WendlandC4{IntType <: Integer, RealType <: AbstractFloat} <: SmoothKernel
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
        kernel::WendlandC4
    )::RealType where RealType <: AbstractFloat

The `kernelValue()` function is used to calculate the kernel value at the distance `r` with the `WendlandC4` kernel.
    
```math
\begin{aligned}
    W(r, h) = \frac{\sigma}{768}
    \begin{cases}
        \begin{array}{ll}
            (2-q)^6(35q^2+36q+12), & 0 \leq q < 2, \\
            0, & q \geq 2,
        \end{array}
    \end{cases}
    \quad \text{where} \quad q = \frac{r}{h}.
\end{aligned}
```
"""
function kernelValue(r::RealType, kernel::WendlandC4)::RealType where {RealType <: AbstractFloat}
    q::RealType = r / kernel.h_
    if q < 2.0
        return kernel.sigma_ * (2.0 - q)^6 * (35.0 * q^2 + 36.0 * q + 12.0) / 768.0
    else
        return 0.0
    end
end

function kernelGradient(r::RealType, kernel::WendlandC4)::RealType where {RealType <: AbstractFloat}
    q::RealType = r / kernel.h_
    if q < 2.0
        return -kernel.sigma_ / kernel.h_ * (35.0 / 96.0 * q^2 + 7.0 / 48.0 * q) * (2.0 - q)^5
    else
        return 0.0
    end
end
