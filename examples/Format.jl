#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/15 13:48:56
  @ license: MIT
  @ description:
 =#

using JuliaFormatter

for path_here in ["CollapseDry", "Cruchaga2D/", "Cruchaga3D/", "NaturalConvectionCavityNormalization/"]
    JuliaFormatter.format(
        joinpath(@__DIR__, path_here),
        indent = 4,
        margin = 150,
        always_for_in = true,
        whitespace_typedefs = true,
        whitespace_ops_in_indices = true,
        remove_extra_newlines = true,
        pipe_to_function_call = false,
        always_use_return = true,
        whitespace_in_kwargs = true,
        trailing_comma = true,
    )
end
