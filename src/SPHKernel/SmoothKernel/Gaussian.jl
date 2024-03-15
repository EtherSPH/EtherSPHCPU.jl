#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/09 17:15:41
  @ license: MIT
  @ description:
 =#

struct Gaussian{IntType <: Integer, RealType <: AbstractFloat} <: SmoothKernel
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
        kernel::Gaussian
    )::RealType where RealType <: AbstractFloat

The `kernelValue()` function is used to calculate the kernel value at the distance `r` with the `Gaussian` kernel.

```math
\begin{aligned}
    W(r, h) = \sigma \exp(-q^2) \quad \text{where} \quad q = \frac{r}{h}.
\end{aligned}
```
"""
function kernelValue(r::RealType, kernel::Gaussian)::RealType where {RealType <: AbstractFloat}
    q::RealType = r / kernel.h_
    if q < 3.0
        return kernel.sigma_ * exp(-q^2)
    else
        return 0.0
    end
end

"""
    kernelGradient(
        r::RealType where RealType <: AbstractFloat,
        kernel::Gaussian
    )::RealType where RealType <: AbstractFloat

The `kernelGradient()` function is used to calculate the kernel gradient at the distance `r` with the `Gaussian` kernel.
See also [kernelGradient](@ref).
"""
function kernelGradient(r::RealType, kernel::Gaussian)::RealType where {RealType <: AbstractFloat}
    q::RealType = r / kernel.h_
    if q < 3.0
        return -2.0 * kernel.sigma_ / kernel.h_ * q * exp(-q^2)
    else
        return 0.0
    end
end
