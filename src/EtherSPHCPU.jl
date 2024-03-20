#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/09 14:27:53
  @ license: MIT
  @ description:
 =#

module EtherSPHCPU # EtherSPHCPU module begin

using LinearAlgebra
using FLoops
FLoops.assistant(:ignore)
using CellListMap
using WriteVTK
using Dates
using ProgressBars
using BangBang
using MicroCollections
using TypeTree
using ExportAll

include("./SPHKernel/SPHKernel.jl")
include("./Particle/Particle.jl")
include("./EquationModel/EquationModel.jl")
include("./Neighbour/Neighbour.jl")
include("./ParticleAction/ParticleAction.jl")
include("./TimeDiscretization/TimeDiscretization.jl")
include("./DataIO/DataIO.jl")
include("./Solve/Solve.jl")

""" 
    showAbstractTypeTree(abstract_type::DataType)::Nothing

The `showAbstractTypeTree()` function is used to show the type tree of the abstract type.
"""
function showAbstractTypeTree(abstract_type::DataType)::Nothing
    for line::String in TypeTree.tt(abstract_type)
        print(line)
    end
    return nothing
end

"""
    __init__()

The `__init__()` function is the first function called when the module is loaded.
It's used to print all the abstract type tree of the module and a welcome message.
"""
function __init__()::Nothing
    println("="^100)
    for ether_sph_julia_cpu_type in
        [AbstractSPHKernel, AbstractParticle, AbstractEquationModel, AbstractNeighbour, AbstractTimeDiscretization, AbstractDataIO]
        showAbstractTypeTree(ether_sph_julia_cpu_type)
        println("-"^80)
    end
    println("Welcome to EtherSPHCPU.jl !")
    println("="^100)
    return nothing
end

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

@exportAll

end # module EtherSPHCPU end
