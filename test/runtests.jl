using Test
using EtherSPHCPU
using JuliaFormatter

JuliaFormatter.format(
    ".",
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

@testset "EtherSPHCPU" begin
    include("./SPHKernel.jl")
    include("./Particle.jl")
    include("./EquationModel.jl")
    include("./Neighbour.jl")
end
