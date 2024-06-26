#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/20 10:13:11
  @ license: MIT
  @ description:
 =#

using EtherSPHCPU

const reynolds_number = 100.0;

const lid_length = 1.0;
const lid_velocity = 1.0;
const rho_0 = 1.0;
const mu_0 = lid_length * rho_0 * lid_velocity / reynolds_number;
const c_0 = 10 * lid_velocity;
const p_0 = 0.001 * rho_0 * c_0^2;
const gamma = 7;

const dim = 2;
const dr = 0.01;
const gap = dr;
const influence_radius = 3.0 * dr;

const smooth_kernel = SmoothKernel(influence_radius, dim, WendlandC2);

const body_force_vec = [0.0, 0.0];

const wc_lm = CommonWeaklyCompressibleLiquidModel(rho_0, c_0, p_0, gamma, mu_0, body_force_vec);

const dt = 0.1 * smooth_kernel.influence_radius_ / c_0;
const total_time = 15.0;
const total_step = total_time / dt |> round |> Int;
const output_interval = 100;
const density_reinitialized_interval = 100;
const check_interval = 100;
const dr_forward_euler =
    DensityReinitializedFixedTimeIntervalForwardEuler(dt, total_step, output_interval, density_reinitialized_interval, check_interval);

const step_digit = div(total_step, output_interval) |> Int |> string |> length;
const file_name = "lid_driven_cavity";
const file_suffix = ".vtp";
const dir_path = "./examples/LidDrivenCavity/LidDrivenCavityData/";
const fileld_symbols = [:rho_, :p_, :c_];
const field_names = ["Density", "Pressure", "SoundSpeed"];
vtp_io = VTPIO(step_digit, file_name, file_suffix, dir_path, fileld_symbols, field_names);

const x0 = 0.0;
const y0 = 0.0;

function createRectangleParticles(particle_type, dx, x_0, y_0, width_number, height_number)
    particles = [particle_type(2) for _ in 1:(width_number * height_number)]
    for i_row in 1:height_number, j_col in 1:width_number
        index = (i_row - 1) * width_number + j_col
        particles[index].x_vec_ = [x_0 + (j_col - 1) * dx, y_0 + (i_row - 1) * dx] .+ dx / 2
    end
    return particles
end

fluid_row_number = lid_length / dr |> round |> Int;
fluid_col_number = lid_length / dr |> round |> Int;
fluid_particles = createRectangleParticles(CommonLiquidParticle, dr, x0, y0, fluid_col_number, fluid_row_number);
for index in eachindex(fluid_particles)
    fluid_particles[index].rho_ = rho_0
    updatePressure!(fluid_particles[index], wc_lm)
    fluid_particles[index].c_ = c_0
    fluid_particles[index].mass_ = fluid_particles[index].rho_ * dr^2
    fluid_particles[index].gap_ = gap
end

const lid_thick_number = 3;
const lid_length_number = lid_length / dr |> round |> Int;

lid_particles = createRectangleParticles(StaticVelocityWallParticle, dr, x0, y0 + lid_length, lid_length_number, lid_thick_number);
for index in eachindex(lid_particles)
    lid_particles[index].v_vec_ = [lid_velocity, 0.0]
    lid_particles[index].gap_ = gap
    lid_particles[index].normal_vec_ = [0.0, -1.0]
end

left_particles = createRectangleParticles(StaticVelocityWallParticle, dr, x0 - lid_thick_number * dr, y0, lid_thick_number, lid_length_number);
for index in eachindex(left_particles)
    left_particles[index].v_vec_ = [0.0, 0.0]
    left_particles[index].gap_ = gap
    left_particles[index].normal_vec_ = [1.0, 0.0]
end

# right StaticVelocityWallParticle
right_particles = createRectangleParticles(StaticVelocityWallParticle, dr, x0 + lid_length, y0, lid_thick_number, lid_length_number);
for index in eachindex(right_particles)
    right_particles[index].v_vec_ = [0.0, 0.0]
    right_particles[index].gap_ = gap
    right_particles[index].normal_vec_ = [-1.0, 0.0]
end

# bottom StaticVelocityWallParticle
bottom_particles = createRectangleParticles(StaticVelocityWallParticle, dr, x0, y0 - lid_thick_number * dr, lid_length_number, lid_thick_number);
for index in eachindex(bottom_particles)
    bottom_particles[index].v_vec_ = [0.0, 0.0]
    bottom_particles[index].gap_ = gap
    bottom_particles[index].normal_vec_ = [0.0, 1.0]
end

# left bottom
left_bottom_particles = createRectangleParticles(
    StaticVelocityWallParticle,
    dr,
    x0 - lid_thick_number * dr,
    y0 - lid_thick_number * dr,
    lid_thick_number,
    lid_thick_number,
);
for index in eachindex(left_bottom_particles)
    left_bottom_particles[index].v_vec_ = [0.0, 0.0]
    left_bottom_particles[index].gap_ = gap
    left_bottom_particles[index].normal_vec_ = [1.0, 1.0] ./ sqrt(2.0)
end

# right bottom
right_bottom_particles =
    createRectangleParticles(StaticVelocityWallParticle, dr, x0 + lid_length, y0 - lid_thick_number * dr, lid_thick_number, lid_thick_number);
for index in eachindex(right_bottom_particles)
    right_bottom_particles[index].v_vec_ = [0.0, 0.0]
    right_bottom_particles[index].gap_ = gap
    right_bottom_particles[index].normal_vec_ = [-1.0, 1.0] ./ sqrt(2.0)
end

# left top
left_top_particles =
    createRectangleParticles(StaticVelocityWallParticle, dr, x0 - lid_thick_number * dr, y0 + lid_length, lid_thick_number, lid_thick_number);
for index in eachindex(left_top_particles)
    left_top_particles[index].v_vec_ = [lid_velocity, 0.0]
    left_top_particles[index].gap_ = gap
    left_top_particles[index].normal_vec_ = [1.0, -1.0] ./ sqrt(2.0)
end

# right top
right_top_particles = createRectangleParticles(StaticVelocityWallParticle, dr, x0 + lid_length, y0 + lid_length, lid_thick_number, lid_thick_number);
for index in eachindex(right_top_particles)
    right_top_particles[index].v_vec_ = [lid_velocity, 0.0]
    right_top_particles[index].gap_ = gap
    right_top_particles[index].normal_vec_ = [-1.0, -1.0] ./ sqrt(2.0)
end

check(particle::AbstractParticle)::Bool = true;

velocity_particles = vcat(
    lid_particles,
    left_particles,
    right_particles,
    bottom_particles,
    left_bottom_particles,
    right_bottom_particles,
    left_top_particles,
    right_top_particles,
);

main() = @fastmath solve!(fluid_particles, velocity_particles, smooth_kernel, wc_lm, dr_forward_euler, vtp_io)
