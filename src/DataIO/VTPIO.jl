#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/10 22:45:09
  @ license: MIT
  @ description:
 =#

struct VTPIO{IntType <: Integer} <: AbstractDataIO
    step_digit_::IntType
    file_name_::String
    file_suffix_::String
    dir_path_::String
    field_symbols_::Vector{Symbol}
    field_names_::Vector{String}
end

function vtpFileNameAtStep(step::IntType where {IntType <: Integer}, vtp_io::VTPIO)::String
    return joinpath(vtp_io.dir_path_, string(vtp_io.file_name_, string(step, pad = vtp_io.step_digit_), vtp_io.file_suffix_))
end

function getParticleFields!(
    field_symbols::Vector{Symbol},
    particles::ParticleArrayType where {ParticleArrayType <: AbstractVector{<:AbstractParticle}},
    fields,
)::Nothing
    particle_field_symbols = fieldnames(eltype(particles))
    for i_field in eachindex(field_symbols)
        if field_symbols[i_field] in particle_field_symbols
            @floop @simd for i_particle in eachindex(particles)
                @inbounds fields[i_field, i_particle] = getfield(particles[i_particle], field_symbols[i_field])
            end
        else
            @floop @simd for i_particle in eachindex(particles)
                @inbounds fields[i_field, i_particle] = NaN
            end
        end
    end
    return nothing
end

function writeVTP(
    step::IntType,
    t::RealType,
    vtp_io::VTPIO,
    particles_group::ParticlesArrayType where {ParticlesArrayType <: AbstractVector{<:AbstractVector}},
)::Nothing where {IntType <: Integer, RealType <: AbstractFloat}
    dim::IntType = length(particles_group[1][1].x_vec_)
    n_particles_group_vector = [length(particles) for particles in particles_group]
    n_particles = sum(n_particles_group_vector)
    points = zeros(RealType, dim, n_particles)
    velocitys = zeros(RealType, dim, n_particles)
    fields = zeros(RealType, length(vtp_io.field_symbols_), n_particles)
    types = zeros(IntType, n_particles)
    current_point_index::IntType = 1
    current_particles_type_index::IntType = 1
    for particles in particles_group
        current_particles_number = length(particles)
        let current_particles_number = current_particles_number
            let current_point_index = current_point_index
                let current_particles_type_index = current_particles_type_index
                    @views getParticleFields!(
                        vtp_io.field_symbols_,
                        particles,
                        fields[:, current_point_index:(current_point_index + current_particles_number - 1)],
                    )
                    @floop @simd for i_particle in eachindex(particles)
                        @inbounds points[:, current_point_index + i_particle - 1] .= particles[i_particle].x_vec_
                        @inbounds types[current_point_index + i_particle - 1] = current_particles_type_index
                    end
                    if eltype(particles) <: MovableParticle || eltype(particles) <: StaticVelocityWallParticle
                        @floop @simd for i_particle in eachindex(particles)
                            @inbounds velocitys[:, current_point_index + i_particle - 1] .= particles[i_particle].v_vec_
                        end
                    else
                        @floop @simd for i_particle in eachindex(particles)
                            @inbounds velocitys[:, current_point_index + i_particle - 1] .= NaN
                        end
                    end
                end
            end
        end
        current_point_index += current_particles_number
        current_particles_type_index += 1
    end
    vtp_file_name = vtpFileNameAtStep(step, vtp_io)
    cells = [MeshCell(PolyData.Verts(), [i]) for i in 1:n_particles]
    vtk_file = vtk_grid(vtp_file_name, points, cells)
    vtk_file["Step"] = step
    vtk_file["Time"] = t
    vtk_file["WallTime"] = Dates.format(Dates.now(), kWallTimeFormat)
    vtk_file["Velocity"] = velocitys
    vtk_file["Type"] = types
    for i_field in eachindex(vtp_io.field_symbols_)
        @inbounds vtk_file[vtp_io.field_names_[i_field]] = fields[i_field, :]
    end
    vtk_save(vtk_file)
    return nothing
end
