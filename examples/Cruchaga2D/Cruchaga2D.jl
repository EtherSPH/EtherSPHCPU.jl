#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/11 22:11:55
  @ license: MIT
  @ description:
 =#

using EtherSPHCPU

# M.A. Cruchaga

const water_width = 0.114
const water_height = 0.114
const dr = 0.002
const gap = dr
const fluid_row_number = water_height / dr |> round |> Int;
const fluid_col_number = water_width / dr |> round |> Int;
const dim = 2;
const influence_radius = 3.0 * dr;

smooth_kernel = SmoothKernel(influence_radius, dim, WendlandC2);

const box_width = 0.42;
const box_height = 0.44;

const gravity = 9.81;
const g_vec = [0.0, -gravity];
const rho_0 = 1000.0;
const c_0 = 15.0;
const p_0 = 0.0;
const gamma = 7;
const mu_0 = 1e-3;
wc_lm = CommonWeaklyCompressibleLiquidModel(rho_0, c_0, p_0, gamma, mu_0, g_vec);

const dt = 0.1 * smooth_kernel.influence_radius_ / c_0;
const total_time = 1.0;
const total_step = total_time / dt |> round |> Int;
const output_interval = 100;
const density_reinitialized_interval = 5;
const check_interval = 10;
dr_fti_forward_euler =
    DensityReinitializedFixedTimeIntervalForwardEuler(dt, total_step, output_interval, check_interval, density_reinitialized_interval);

const step_digit = div(total_step, output_interval) |> Int |> string |> length;
const file_name = "cruchaga_2d";
const file_suffix = ".vtp";
const dir_path = "./examples/Cruchaga2D/Cruchaga2DData";
const fileld_symbols = [:rho_, :p_, :c_];
const fileld_names = ["Density", "Pressure", "SoundSpeed"];
vtp_io = VTPIO(step_digit, file_name, file_suffix, dir_path, fileld_symbols, fileld_names);

const x0 = 0.0;
const y0 = 0.0;
const fluid_number = fluid_row_number * fluid_col_number;

function createRectangleParticles(partricle_type, dx, x_0, y_0, width_number, height_number)
    particles = [partricle_type(2) for _ in 1:(width_number * height_number)]
    for i_row in 1:height_number, j_col in 1:width_number
        index = (i_row - 1) * width_number + j_col
        particles[index].x_vec_ .= [x_0 + (j_col - 1) * dx, y_0 + (i_row - 1) * dx] .+ dx / 2
    end
    return particles
end

fluid_particles = createRectangleParticles(CommonLiquidParticle, dr, x0, y0, fluid_col_number, fluid_row_number);
for index in eachindex(fluid_particles)
    fluid_particles[index].p_ = rho_0 * gravity * (water_height - fluid_particles[index].x_vec_[2])
    fluid_particles[index].rho_ = (fluid_particles[index].p_ / wc_lm.b_ + 1.0)^(1.0 / gamma) * rho_0
    updatePressure!(fluid_particles[index], wc_lm)
    fluid_particles[index].mass_ = dr^2 * fluid_particles[index].rho_
    fluid_particles[index].gap_ = dr
end

const wall_thick_number = 3;
const wall_box_width_number = box_width / dr |> round |> Int;
const wall_box_height_number = box_height / dr |> round |> Int;

bottom_wall_particles =
    createRectangleParticles(CompulsiveWallParticle, dr, x0, y0 - dr * wall_thick_number, wall_box_width_number, wall_thick_number);
for index in eachindex(bottom_wall_particles)
    bottom_wall_particles[index].normal_vec_ .= [0.0, 1.0]
    bottom_wall_particles[index].gap_ = dr
end

left_wall_particles = createRectangleParticles(CompulsiveWallParticle, dr, x0 - dr * wall_thick_number, y0, wall_thick_number, wall_box_height_number);
for index in eachindex(left_wall_particles)
    left_wall_particles[index].normal_vec_ .= [1.0, 0.0]
    left_wall_particles[index].gap_ = dr
end

right_wall_particles = createRectangleParticles(CompulsiveWallParticle, dr, x0 + box_width, y0, wall_thick_number, wall_box_height_number);
for index in eachindex(right_wall_particles)
    right_wall_particles[index].normal_vec_ .= [-1.0, 0.0]
    right_wall_particles[index].gap_ = dr
end

top_wall_particles = createRectangleParticles(CompulsiveWallParticle, dr, x0, y0 + box_height, wall_box_width_number, wall_thick_number);
for index in eachindex(top_wall_particles)
    top_wall_particles[index].normal_vec_ .= [0.0, -1.0]
    top_wall_particles[index].gap_ = dr
end

left_bottom_corner_particle = createRectangleParticles(
    CompulsiveWallParticle,
    dr,
    x0 - dr * wall_thick_number,
    y0 - dr * wall_thick_number,
    wall_thick_number,
    wall_thick_number,
);
for index in eachindex(left_bottom_corner_particle)
    left_bottom_corner_particle[index].normal_vec_ .= [1.0, 1.0] ./ sqrt(2)
    left_bottom_corner_particle[index].gap_ = dr
end

right_bottom_corner_particle =
    createRectangleParticles(CompulsiveWallParticle, dr, x0 + box_width, y0 - dr * wall_thick_number, wall_thick_number, wall_thick_number);
for index in eachindex(right_bottom_corner_particle)
    right_bottom_corner_particle[index].normal_vec_ .= [-1.0, 1.0] ./ sqrt(2)
    right_bottom_corner_particle[index].gap_ = dr
end

left_top_corner_particle =
    createRectangleParticles(CompulsiveWallParticle, dr, x0 - dr * wall_thick_number, y0 + box_height, wall_thick_number, wall_thick_number);

for index in eachindex(left_top_corner_particle)
    left_top_corner_particle[index].normal_vec_ .= [1.0, -1.0] ./ sqrt(2)
    left_top_corner_particle[index].gap_ = dr
end

right_top_corner_particle =
    createRectangleParticles(CompulsiveWallParticle, dr, x0 + box_width, y0 + box_height, wall_thick_number, wall_thick_number);

for index in eachindex(right_top_corner_particle)
    right_top_corner_particle[index].normal_vec_ .= [-1.0, -1.0] ./ sqrt(2)
    right_top_corner_particle[index].gap_ = dr
end

wall_particles = vcat(
    bottom_wall_particles,
    left_wall_particles,
    right_wall_particles,
    top_wall_particles,
    left_bottom_corner_particle,
    right_bottom_corner_particle,
    left_top_corner_particle,
    right_top_corner_particle,
);

main() = @fastmath solve!(fluid_particles, wall_particles, smooth_kernel, wc_lm, dr_fti_forward_euler, vtp_io);
