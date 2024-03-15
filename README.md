# EtherSPHCPU.jl

Julia parallel version of EtherSPH on CPU.

The main idea of `EtherSPHCPU` is to take advantage of julia's multi-dispatch feature to apply different behaviours to different types of particles, kernels, models, etc. This is achieved by defining abstract types and then defining concrete types that inherit from the abstract types. You may customize your own types to add more equations and physic models with this structure.

## Tests

Run a julia REPL and activate the environment.

```bash
xxx/xxx/EtherSPHCPU.jl> julia
```

After that, you may run the following commands:

```julia
pkg> activate .
pkg> instantiate # if you haven't instantiate the environment
pkg> resolve # if you haven't resolve the environment
pkg> test
```

## Docs

```bash
julia --project=docs docs/make.jl # currently deployed locally
```

## Examples

```bash
pkg> activate examples
julia> include("examples/xxx/xxx.jl")
julia> @time main()
```

You may also run `Pluto` to monitor the simulation progress or do post-processing as each example is acompanied with a `xxxPost.jl` file.

## Data Type Sturcture

Thanks to [TypeTree.jl](https://github.com/cnaak/TypeTree.jl), as `EtherSPHCPU` is introduced, the data type structure is shown as follows:

```bash
====================================================================================================
EtherSPHCPU.AbstractSPHKernel
 └─ EtherSPHCPU.SmoothKernel
     ├─ EtherSPHCPU.CubicSpline
     ├─ EtherSPHCPU.Gaussian
     ├─ EtherSPHCPU.WendlandC2
     └─ EtherSPHCPU.WendlandC4
--------------------------------------------------------------------------------
EtherSPHCPU.AbstractParticle
 ├─ EtherSPHCPU.FixedParticle
 │   └─ EtherSPHCPU.WallParticle
 │       ├─ EtherSPHCPU.CompulsiveWallParticle
 │       ├─ EtherSPHCPU.StaticVelocityWallParticle
 │       └─ EtherSPHCPU.ThermostaticWallParticle
 └─ EtherSPHCPU.MovableParticle
     └─ EtherSPHCPU.FluidParticle
         └─ EtherSPHCPU.LiquidParticle
             ├─ EtherSPHCPU.CommonLiquidParticle
             └─ EtherSPHCPU.ThermalLiquidParticle
--------------------------------------------------------------------------------
EtherSPHCPU.AbstractEquationModel
 └─ EtherSPHCPU.LiquidModel
     └─ EtherSPHCPU.WeaklyCompressibleLiquidModel
         ├─ EtherSPHCPU.CommonWeaklyCompressibleLiquidModel
         └─ EtherSPHCPU.ThermalWeaklyCompressibleLiquidModel
--------------------------------------------------------------------------------
EtherSPHCPU.AbstractNeighbour
 └─ EtherSPHCPU.CommonNeighbour
--------------------------------------------------------------------------------
EtherSPHCPU.AbstractTimeDiscretization
 ├─ EtherSPHCPU.FixedTimeIntervalDiscretization
 │   └─ EtherSPHCPU.DensityReinitializedFixedTimeIntervalForwardEuler
 └─ EtherSPHCPU.VariableTimeIntervalDiscretization
--------------------------------------------------------------------------------
EtherSPHCPU.AbstractDataIO
 └─ EtherSPHCPU.VTPIO
--------------------------------------------------------------------------------
Welcome to EtherSPHCPU.jl !
====================================================================================================
```

## Thanks

- [TypeTree.jl](https://github.com/cnaak/TypeTree.jl): Never let me get lost in my abstract types.
- [CellListMap.jl](https://github.com/m3g/CellListMap.jl): The core of the neighbour search algorithm. It seems the author is a researcher on molecular dynamics.
- [FLoops.jl](git@github.com:JuliaFolds/FLoops.jl.git): Fast parallel cumputing with a single marco `@`! It speeds up around 2-3 times compared with `Threads.@threads` during my tests.
- [SmoothedParticles.jl](https://github.com/OndrejKincl/SmoothedParticles.jl): This package inspires me at the beginning of the project. It has good examples and a good structure. I learned a lot from it.