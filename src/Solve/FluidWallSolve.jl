#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/11 21:19:03
  @ license: MIT
  @ description:
 =#

function eachStep!(
    fluid_particles::FluidParticlesType where {FluidParticlesType <: AbstractVector{<:FluidParticle}},
    wall_particles::WallParticlesType where {WallParticlesType <: AbstractVector{<:WallParticle}},
    smooth_kernel::SmoothKernelType where {SmoothKernelType <: SmoothKernel},
    wc_lm::WeaklyCompressibleLiquidModelType where {WeaklyCompressibleLiquidModelType <: WeaklyCompressibleLiquidModel},
    dr_fti_forward_euler::DensityReinitializedFixedTimeIntervalForwardEuler,
    fluid_neighbour_system::InPlaceNeighborList,
    fluid_wall_neighbour_system::InPlaceNeighborList,
    fluid_neighbours::NeighbourListType where {NeighbourListType <: AbstractVector{<:AbstractNeighbour}},
    fluid_wall_neighbours::NeighbourListType where {NeighbourListType <: AbstractVector{<:AbstractNeighbour}},
    step::IntType where {IntType <: Integer};
    check::Function = defaultCheck,
)::Nothing
    fluid_neighbours_task = @async begin
        update!(fluid_neighbour_system, [p.x_vec_ for p in fluid_particles]; cutoff = smooth_kernel.influence_radius_)
        findNeighbours!(neighborlist!(fluid_neighbour_system), fluid_neighbours, fluid_particles, fluid_particles, smooth_kernel)
    end
    fluid_wall_neighbours_task = @async begin
        update!(
            fluid_wall_neighbour_system,
            [p.x_vec_ for p in fluid_particles],
            [p.x_vec_ for p in wall_particles];
            cutoff = smooth_kernel.influence_radius_,
        )
        findNeighbours!(neighborlist!(fluid_wall_neighbour_system), fluid_wall_neighbours, fluid_particles, wall_particles, smooth_kernel)
    end
    wait(fluid_neighbours_task)
    wait(fluid_wall_neighbours_task)
    # continuity equation
    @floop @simd for fluid_neighbour in fluid_neighbours
        @inbounds continuity!(fluid_particles[fluid_neighbour.i_], fluid_particles[fluid_neighbour.j_], fluid_neighbour)
    end
    # update density and pressure
    @floop @simd for fluid_particle in fluid_particles
        updateDensity!(fluid_particle, dr_fti_forward_euler.dt_)
        updatePressure!(fluid_particle, wc_lm)
    end
    # momentum equation
    fluid_force_task = @async begin
        @floop @simd for fluid_neighbour in fluid_neighbours
            @inbounds pressureForce!(fluid_particles[fluid_neighbour.i_], fluid_particles[fluid_neighbour.j_], fluid_neighbour, smooth_kernel)
            @inbounds viscosityForce!(fluid_particles[fluid_neighbour.i_], fluid_particles[fluid_neighbour.j_], fluid_neighbour, smooth_kernel, wc_lm)
        end
    end
    fluid_wall_force_task = @async begin
        @floop @simd for fluid_wall_neighbour in fluid_wall_neighbours
            @inbounds wallForce!(fluid_particles[fluid_wall_neighbour.i_], wall_particles[fluid_wall_neighbour.j_], fluid_wall_neighbour, smooth_kernel)
            @inbounds viscosityForce!(
                fluid_particles[fluid_wall_neighbour.i_],
                wall_particles[fluid_wall_neighbour.j_],
                fluid_wall_neighbour,
                smooth_kernel,
                wc_lm,
            )
        end
    end
    wait(fluid_force_task)
    wait(fluid_wall_force_task)
    @floop @simd for fluid_particle in fluid_particles
        updateVelocity!(fluid_particle, dr_fti_forward_euler.dt_, wc_lm.body_force_vec_)
        updatePosition!(fluid_particle, dr_fti_forward_euler.dt_)
    end
    if isDensityReinitializedStep(step, dr_fti_forward_euler)
        reconstructScalar!(fluid_particles, :rho_, fluid_neighbours, smooth_kernel)
    end
    if isCheckStep(step, dr_fti_forward_euler)
        @floop @simd for i_particle in eachindex(fluid_particles)
            out_of_bounds_here = SingletonVector(i_particle)
            if !check(fluid_particles[i_particle])
                @reduce(out_of_bounds_global = append!!(EmptyVector(), out_of_bounds_here))
            end
        end
        deleteat!(fluid_particles, out_of_bounds_global)
    end
    return nothing
end

function solve!(
    fluid_particles::FluidParticlesType where {FluidParticlesType <: AbstractVector{<:FluidParticle}},
    wall_particles::WallParticlesType where {WallParticlesType <: AbstractVector{<:WallParticle}},
    smooth_kernel::SmoothKernelType where {SmoothKernelType <: SmoothKernel},
    wc_lm::WeaklyCompressibleLiquidModelType where {WeaklyCompressibleLiquidModelType <: WeaklyCompressibleLiquidModel},
    dr_fti_forward_euler::DensityReinitializedFixedTimeIntervalForwardEuler,
    vtp_io::VTPIO;
    check::Function = defaultCheck,
)::Nothing
    assureDirPathExist(vtp_io)
    writeVTP(0, 0.0, vtp_io, [fluid_particles, wall_particles])
    fluid_neighbour_system = InPlaceNeighborList(x = [p.x_vec_ for p in fluid_particles], cutoff = smooth_kernel.influence_radius_, parallel = true)
    fluid_wall_neighbour_system = InPlaceNeighborList(
        x = [p.x_vec_ for p in fluid_particles],
        y = [p.x_vec_ for p in wall_particles],
        cutoff = smooth_kernel.influence_radius_,
        parallel = true,
    )
    fluid_neighbours = CommonNeighbour[] # to be initialized
    fluid_wall_neighbours = CommonNeighbour[] # to be initialized
    for step in ProgressBar(1:(dr_fti_forward_euler.total_step_))
        # eachStep!(
        #     fluid_particles, wall_particles,
        #     smooth_kernel,
        #     wc_lm,
        #     dr_fti_forward_euler,
        #     fluid_neighbour_system, fluid_wall_neighbour_system,
        #     fluid_neighbours, fluid_wall_neighbours,
        #     step;
        #     check=check
        # )
        try
            eachStep!(
                fluid_particles,
                wall_particles,
                smooth_kernel,
                wc_lm,
                dr_fti_forward_euler,
                fluid_neighbour_system,
                fluid_wall_neighbour_system,
                fluid_neighbours,
                fluid_wall_neighbours,
                step;
                check = check,
            )
        catch e
            println(e)
            println("Error at step: ", step)
            writeVTP(step, step * dr_fti_forward_euler.dt_, vtp_io, [fluid_particles, wall_particles])
        end
        if isOutputStep(step, dr_fti_forward_euler)
            writeVTP(div(step, dr_fti_forward_euler.output_interval_), step * dr_fti_forward_euler.dt_, vtp_io, [fluid_particles, wall_particles])
        end
    end
    return nothing
end
