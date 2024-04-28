#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/04/25 21:05:52
  @ license: MIT
  @ description:
 =#

function eachStep!(
    fluid_particles::FluidParticlesType where {FluidParticlesType <: AbstractVector{<:FluidParticle}},
    inlet_fluid_particles::FluidParticlesType where {FluidParticlesType <: AbstractVector{<:FluidParticle}},
    outlet_fluid_particles::FluidParticlesType where {FluidParticlesType <: AbstractVector{<:FluidParticle}},
    wall_particles::WallParticlesType where {WallParticlesType <: AbstractVector{<:WallParticle}},
    smooth_kernel::SmoothKernelType where {SmoothKernelType <: SmoothKernel},
    wc_lm::WeaklyCompressibleLiquidModelType where {WeaklyCompressibleLiquidModelType <: WeaklyCompressibleLiquidModel},
    dr_fti_forward_euler::DensityReinitializedFixedTimeIntervalForwardEuler,
    fluid_neighbour_system::InPlaceNeighborList,
    fluid_wall_neighbour_system::InPlaceNeighborList,
    fluid_inlet_neighbour_system::InPlaceNeighborList,
    fluid_outlet_neighbour_system::InPlaceNeighborList,
    fluid_neighbours::NeighbourListType where {NeighbourListType <: AbstractVector{<:AbstractNeighbour}},
    fluid_wall_neighbours::NeighbourListType where {NeighbourListType <: AbstractVector{<:AbstractNeighbour}},
    fluid_inlet_neighbours::NeighbourListType where {NeighbourListType <: AbstractVector{<:AbstractNeighbour}},
    fluid_outlet_neighbours::NeighbourListType where {NeighbourListType <: AbstractVector{<:AbstractNeighbour}},
    step::IntType where {IntType <: Integer};
    inlet_into_calculation_domain_check::Function = defaultCheck,
    outlet_out_calculation_domain_check::Function = defaultCheck,
    outlet_out_check::Function = defaultCheck,
    inlet_into_calculation_domain_modify::Function = defaultModify,
    user_defined_function::Function = defaultUserDefinedFunction,
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
    fluid_inlet_neighbours_task = @async begin
        update!(
            fluid_inlet_neighbour_system,
            [p.x_vec_ for p in fluid_particles],
            [p.x_vec_ for p in inlet_fluid_particles];
            cutoff = smooth_kernel.influence_radius_,
        )
        findNeighbours!(neighborlist!(fluid_inlet_neighbour_system), fluid_inlet_neighbours, fluid_particles, inlet_fluid_particles, smooth_kernel)
    end
    fluid_outlet_neighbours_task = @async begin
        update!(
            fluid_outlet_neighbour_system,
            [p.x_vec_ for p in fluid_particles],
            [p.x_vec_ for p in outlet_fluid_particles];
            cutoff = smooth_kernel.influence_radius_,
        )
        findNeighbours!(neighborlist!(fluid_outlet_neighbour_system), fluid_outlet_neighbours, fluid_particles, outlet_fluid_particles, smooth_kernel)
    end
    wait(fluid_neighbours_task)
    wait(fluid_wall_neighbours_task)
    wait(fluid_inlet_neighbours_task)
    wait(fluid_outlet_neighbours_task)
    # continuity equation
    fluid_continuity_task = @async begin
        @floop @simd for fluid_neighbour in fluid_neighbours
            @inbounds continuity!(fluid_particles[fluid_neighbour.i_], fluid_particles[fluid_neighbour.j_], fluid_neighbour)
        end
    end
    inlet_continuity_task = @async begin
        @floop @simd for fluid_inlet_neighbour in fluid_inlet_neighbours
            @inbounds continuity!(fluid_particles[fluid_inlet_neighbour.i_], inlet_fluid_particles[fluid_inlet_neighbour.j_], fluid_inlet_neighbour)
        end
    end
    outlet_continuity_task = @async begin
        @floop @simd for fluid_outlet_neighbour in fluid_outlet_neighbours
            @inbounds continuity!(fluid_particles[fluid_outlet_neighbour.i_], outlet_fluid_particles[fluid_outlet_neighbour.j_], fluid_outlet_neighbour)
        end
    end
    wait(fluid_continuity_task)
    wait(inlet_continuity_task)
    wait(outlet_continuity_task)
    # update density and pressure
    fluid_density_pressure_task = @async begin
        @floop @simd for fluid_particle in fluid_particles
            updateDensity!(fluid_particle, dr_fti_forward_euler.dt_)
            updatePressure!(fluid_particle, wc_lm)
        end
    end
    inlet_density_pressure_task = @async begin
        @floop @simd for inlet_fluid_particle in inlet_fluid_particles
            inlet_fluid_particle.drho_[] = 0.0
        end
    end
    outlet_density_pressure_task = @async begin
        @floop @simd for outlet_fluid_particle in outlet_fluid_particles
            outlet_fluid_particle.drho_[] = 0.0
        end
    end
    wait(fluid_density_pressure_task)
    wait(inlet_density_pressure_task)
    wait(outlet_density_pressure_task)
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
            @inbounds pressureForce!(fluid_particles[fluid_wall_neighbour.i_], fluid_wall_neighbour, smooth_kernel, wc_lm)
        end
    end
    fluid_inlet_force_task = @async begin
        @floop @simd for fluid_inlet_neighbour in fluid_inlet_neighbours
            @inbounds pressureForce!(
                fluid_particles[fluid_inlet_neighbour.i_],
                inlet_fluid_particles[fluid_inlet_neighbour.j_],
                fluid_inlet_neighbour,
                smooth_kernel,
            )
            @inbounds viscosityForce!(
                fluid_particles[fluid_inlet_neighbour.i_],
                inlet_fluid_particles[fluid_inlet_neighbour.j_],
                fluid_inlet_neighbour,
                smooth_kernel,
                wc_lm,
            )
        end
    end
    fluid_outlet_force_task = @async begin
        @floop @simd for fluid_outlet_neighbour in fluid_outlet_neighbours
            @inbounds pressureForce!(
                fluid_particles[fluid_outlet_neighbour.i_],
                outlet_fluid_particles[fluid_outlet_neighbour.j_],
                fluid_outlet_neighbour,
                smooth_kernel,
            )
            @inbounds viscosityForce!(
                fluid_particles[fluid_outlet_neighbour.i_],
                outlet_fluid_particles[fluid_outlet_neighbour.j_],
                fluid_outlet_neighbour,
                smooth_kernel,
                wc_lm,
            )
        end
    end
    wait(fluid_force_task)
    wait(fluid_wall_force_task)
    wait(fluid_inlet_force_task)
    wait(fluid_outlet_force_task)
    # update velocity and position
    fluid_velocity_position_task = @async begin
        @floop @simd for fluid_particle in fluid_particles
            updateVelocity!(fluid_particle, dr_fti_forward_euler.dt_, wc_lm.body_force_vec_)
            updatePosition!(fluid_particle, dr_fti_forward_euler.dt_)
        end
    end
    inlet_velocity_position_task = @async begin
        @floop @simd for inlet_fluid_particle in inlet_fluid_particles
            @simd for i_dim in eachindex(inlet_fluid_particle.x_vec_)
                inlet_fluid_particle.dv_vec_[i_dim][] = 0.0
            end
            updatePosition!(inlet_fluid_particle, dr_fti_forward_euler.dt_)
        end
    end
    outlet_velocity_position_task = @async begin
        @floop @simd for outlet_fluid_particle in outlet_fluid_particles
            @simd for i_dim in eachindex(outlet_fluid_particle.x_vec_)
                outlet_fluid_particle.dv_vec_[i_dim][] = 0.0
            end
            updatePosition!(outlet_fluid_particle, dr_fti_forward_euler.dt_)
        end
    end
    wait(fluid_velocity_position_task)
    wait(inlet_velocity_position_task)
    wait(outlet_velocity_position_task)
    if isDensityReinitializedStep(step, dr_fti_forward_euler)
        reconstructScalar!(fluid_particles, :rho_, fluid_neighbours, smooth_kernel)
    end
    # process for inlet / outlet
    # * 1. get the particles that are into / out from inlet / calculation domain / outlet
    # * - inlet into calculation domain's particle ids
    # * - calculation domain into outlet's particle ids
    # * - outlet out outlet particle ids
    inlet_into_calculation_domain_check_task = @async begin
        @floop @simd for i_inlet_particle in eachindex(inlet_fluid_particles)
            inlet_into_calculation_domain_local_ids = SingletonVector(i_inlet_particle)
            if inlet_into_calculation_domain_check(inlet_fluid_particles[i_inlet_particle]) == true
                @reduce(inlet_into_calculation_domain_global_ids = append!!(EmptyVector(), inlet_into_calculation_domain_local_ids))
            end
        end
        inlet_into_calculation_domain_global_ids
    end
    outlet_out_calculation_domain_check_task = @async begin
        @floop @simd for i_fluid_particle in eachindex(fluid_particles)
            outlet_out_calculation_domain_local_ids = SingletonVector(i_fluid_particle)
            if outlet_out_calculation_domain_check(fluid_particles[i_fluid_particle]) == true
                @reduce(outlet_out_calculation_domain_global_ids = append!!(EmptyVector(), outlet_out_calculation_domain_local_ids))
            end
        end
        outlet_out_calculation_domain_global_ids
    end
    outlet_out_check_task = @async begin
        @floop @simd for i_outlet_particle in eachindex(outlet_fluid_particles)
            outlet_out_local_ids = SingletonVector(i_outlet_particle)
            if outlet_out_check(outlet_fluid_particles[i_outlet_particle]) == true
                @reduce(outlet_out_global_ids = append!!(EmptyVector(), outlet_out_local_ids))
            end
        end
        outlet_out_global_ids
    end
    # use fetch to get the global ids instead of wait
    inlet_into_calculation_domain_global_ids = fetch(inlet_into_calculation_domain_check_task)
    outlet_out_calculation_domain_global_ids = fetch(outlet_out_calculation_domain_check_task)
    outlet_out_global_ids = fetch(outlet_out_check_task)
    # * 2. copy and push
    # * - inlet into calculation domain
    # * - calculation domain into outlet
    inlet_into_calculation_domain_copy_and_push_task = @async begin
        for i_inlet_particle in inlet_into_calculation_domain_global_ids
            push!(fluid_particles, deepcopy(inlet_fluid_particles[i_inlet_particle]))
        end
    end
    calculation_domain_into_outlet_copy_and_push_task = @async begin
        for i_fluid_particle in outlet_out_calculation_domain_global_ids
            push!(outlet_fluid_particles, deepcopy(fluid_particles[i_fluid_particle]))
        end
    end
    wait(inlet_into_calculation_domain_copy_and_push_task)
    wait(calculation_domain_into_outlet_copy_and_push_task)
    # * 3. delete calculation domain / outlet particles and modify the inlet particles (reset)
    # * - remove calculation domain particles into outlet
    # * - remove outlet particles out from outlet
    # * - modify inlet particles
    out_of_calculation_domain_delete_task = @async begin
        deleteat!(fluid_particles, outlet_out_calculation_domain_global_ids)
    end
    out_of_outlet_delete_task = @async begin
        deleteat!(outlet_fluid_particles, outlet_out_global_ids)
    end
    modify_inlet_task = @async begin
        @floop @simd for i_inlet_particle in inlet_into_calculation_domain_global_ids
            inlet_into_calculation_domain_modify(inlet_fluid_particles[i_inlet_particle])
        end
    end
    wait(out_of_calculation_domain_delete_task)
    wait(out_of_outlet_delete_task)
    wait(modify_inlet_task)
    user_defined_function(; t = step * dr_fti_forward_euler.dt_)
    return nothing
end

function solve!(
    fluid_particles::FluidParticlesType where {FluidParticlesType <: AbstractVector{<:FluidParticle}},
    inlet_fluid_particles::FluidParticlesType where {FluidParticlesType <: AbstractVector{<:FluidParticle}},
    outlet_fluid_particles::FluidParticlesType where {FluidParticlesType <: AbstractVector{<:FluidParticle}},
    wall_particles::WallParticlesType where {WallParticlesType <: AbstractVector{<:WallParticle}},
    smooth_kernel::SmoothKernelType where {SmoothKernelType <: SmoothKernel},
    wc_lm::WeaklyCompressibleLiquidModelType where {WeaklyCompressibleLiquidModelType <: WeaklyCompressibleLiquidModel},
    dr_fti_forward_euler::DensityReinitializedFixedTimeIntervalForwardEuler,
    vtp_io::VTPIO;
    inlet_into_calculation_domain_check::Function = defaultCheck,
    outlet_out_calculation_domain_check::Function = defaultCheck,
    outlet_out_check::Function = defaultCheck,
    inlet_into_calculation_domain_modify::Function = defaultModify,
    user_defined_function::Function = defaultUserDefinedFunction,
)::Nothing
    assureDirPathExist(vtp_io)
    writeVTP(0, 0.0, vtp_io, [fluid_particles, inlet_fluid_particles, outlet_fluid_particles, wall_particles])
    fluid_neighbour_system = InPlaceNeighborList(x = [p.x_vec_ for p in fluid_particles], cutoff = smooth_kernel.influence_radius_, parallel = true)
    fluid_wall_neighbour_system = InPlaceNeighborList(
        x = [p.x_vec_ for p in fluid_particles],
        y = [p.x_vec_ for p in wall_particles],
        cutoff = smooth_kernel.influence_radius_,
        parallel = true,
    )
    fluid_inlet_neighbour_system = InPlaceNeighborList(
        x = [p.x_vec_ for p in fluid_particles],
        y = [p.x_vec_ for p in inlet_fluid_particles],
        cutoff = smooth_kernel.influence_radius_,
        parallel = true,
    )
    fluid_outlet_neighbour_system = InPlaceNeighborList(
        x = [p.x_vec_ for p in fluid_particles],
        y = [p.x_vec_ for p in outlet_fluid_particles],
        cutoff = smooth_kernel.influence_radius_,
        parallel = true,
    )
    fluid_neighbours = CommonNeighbour[] # to be initialized
    fluid_wall_neighbours = CommonNeighbour[] # to be initialized
    fluid_inlet_neighbours = CommonNeighbour[] # to be initialized
    fluid_outlet_neighbours = CommonNeighbour[] # to be initialized
    for step in ProgressBar(1:(dr_fti_forward_euler.total_step_))
        eachStep!(
            fluid_particles,
            inlet_fluid_particles,
            outlet_fluid_particles,
            wall_particles,
            smooth_kernel,
            wc_lm,
            dr_fti_forward_euler,
            fluid_neighbour_system,
            fluid_wall_neighbour_system,
            fluid_inlet_neighbour_system,
            fluid_outlet_neighbour_system,
            fluid_neighbours,
            fluid_wall_neighbours,
            fluid_inlet_neighbours,
            fluid_outlet_neighbours,
            step;
            inlet_into_calculation_domain_check = inlet_into_calculation_domain_check,
            outlet_out_calculation_domain_check = outlet_out_calculation_domain_check,
            outlet_out_check = outlet_out_check,
            inlet_into_calculation_domain_modify = inlet_into_calculation_domain_modify,
            user_defined_function = user_defined_function,
        )
        if isOutputStep(step, dr_fti_forward_euler)
            writeVTP(
                div(step, dr_fti_forward_euler.output_interval_),
                step * dr_fti_forward_euler.dt_,
                vtp_io,
                [fluid_particles, inlet_fluid_particles, outlet_fluid_particles, wall_particles],
            )
        end
    end
    return nothing
end
