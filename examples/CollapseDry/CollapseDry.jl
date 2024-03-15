#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/13 12:59:40
  @ license: MIT
  @ description:
 =#

using EtherSPHCPU

const dim = 2;
const water_width = 1.0;
const water_height = 2.0;
const dr = 1e-2;
const gap = dr;
const influence_radius = 3.0 * dr;

smooth_kernel = SmoothKernel(influence_radius, dim, CubicSpline);

smooth_kernel = SmoothKernel(influence_radius, dim, CubicSpline);

const box_width = 4.0;
const box_height = 3.0;

const gravity = 9.8;
const g_vec = [0.0, -gravity];
const rho_0 = 1000.0;
const c_0 = 120.0;
const p_0 = 0.0;
const gamma = 7;
const mu_0 = 1e-3;
wc_lm = CommonWeaklyCompressibleLiquidModel(rho_0, c_0, p_0, gamma, mu_0, g_vec);

const dt = 5e-5;
const total_time = 3.0;
const total_step = total_time / dt |> round |> Int;
const output_interval = 100;
const density_reinitialized_interval = 5
const check_interval = 5;
dr_fti_forward_euler =
    DensityReinitializedFixedTimeIntervalForwardEuler(dt, total_step, output_interval, density_reinitialized_interval, check_interval);

const step_digit = div(total_step, output_interval) |> Int |> string |> length;
const file_name = "collapse_dry";
const file_suffix = ".vtp";
const dir_path = "./examples/CollapseDry/CollapseDryData/";
const fileld_symbols = [:rho_, :p_, :c_];
const fileld_names = ["Density", "Pressure", "SoundSpeed"];
vtp_io = VTPIO(step_digit, file_name, file_suffix, dir_path, fileld_symbols, fileld_names);

const x0 = 0.0;
const y0 = 0.0;

function createRectangleParticles(partricle_type, dx, x_0, y_0, width_number, height_number)
    particles = [partricle_type(2) for _ in 1:(width_number * height_number)]
    for i_row in 1:height_number, j_col in 1:width_number
        index = (i_row - 1) * width_number + j_col
        particles[index].x_vec_ .= [x_0 + (j_col - 1) * dx, y_0 + (i_row - 1) * dx] .+ dx / 2
    end
    return particles
end

const fluid_row_number = water_height / dr |> round |> Int;
const fluid_col_number = water_width / dr |> round |> Int;
const fluid_number = fluid_row_number * fluid_col_number;

fluid_particles = createRectangleParticles(CommonLiquidParticle, dr, x0, y0, fluid_col_number, fluid_row_number);
for index in 1:fluid_number
    fluid_particles[index].p_ = rho_0 * gravity * (water_height + dr - fluid_particles[index].x_vec_[2])
    fluid_particles[index].rho_ = (fluid_particles[index].p_ / wc_lm.b_ + 1.0)^(1.0 / gamma) * rho_0
    # fluid_particles[index].rho_ = rho_0;
    updatePressure!(fluid_particles[index], wc_lm)
    fluid_particles[index].mass_ = fluid_particles[index].rho_ * dr^2
    fluid_particles[index].gap_ = dr
end

const wall_thick_number = 3;
const wall_box_width_number = box_width / dr |> round |> Int;
const wall_box_height_number = box_height / dr |> round |> Int;

bottom_wall_particles =
    createRectangleParticles(CompulsiveWallParticle, dr, x0, y0 - wall_thick_number * dr, wall_box_width_number, wall_thick_number);
for index in eachindex(bottom_wall_particles)
    bottom_wall_particles[index].normal_vec_ = [0.0, 1.0]
    bottom_wall_particles[index].gap_ = dr
end

left_wall_particles = createRectangleParticles(CompulsiveWallParticle, dr, x0 - wall_thick_number * dr, y0, wall_thick_number, wall_box_height_number);
for index in eachindex(left_wall_particles)
    left_wall_particles[index].normal_vec_ = [1.0, 0.0]
    left_wall_particles[index].gap_ = dr
end

right_wall_particles = createRectangleParticles(CompulsiveWallParticle, dr, x0 + box_width, y0, wall_thick_number, wall_box_height_number);
for index in eachindex(right_wall_particles)
    right_wall_particles[index].normal_vec_ = [-1.0, 0.0]
    right_wall_particles[index].gap_ = dr
end

left_bottom_corner_particles = createRectangleParticles(
    CompulsiveWallParticle,
    dr,
    x0 - wall_thick_number * dr,
    y0 - wall_thick_number * dr,
    wall_thick_number,
    wall_thick_number,
);
for index in eachindex(left_bottom_corner_particles)
    left_bottom_corner_particles[index].normal_vec_ = [1.0, 1.0] ./ sqrt(2.0)
    left_bottom_corner_particles[index].gap_ = dr
end

right_bottom_corner_particles =
    createRectangleParticles(CompulsiveWallParticle, dr, x0 + box_width, y0 - wall_thick_number * dr, wall_thick_number, wall_thick_number);
for index in eachindex(right_bottom_corner_particles)
    right_bottom_corner_particles[index].normal_vec_ = [-1.0, 1.0] ./ sqrt(2.0)
    right_bottom_corner_particles[index].gap_ = dr
end

wall_particles = vcat(bottom_wall_particles, left_wall_particles, right_wall_particles, left_bottom_corner_particles, right_bottom_corner_particles);

function check(particle::LiquidParticle)::Bool
    if particle.x_vec_[1] > x0 + box_width
        return false
    elseif particle.x_vec_[2] > y0 + box_height * 1.0
        return false
    else
        return true
    end
end

main() = @fastmath solve!(fluid_particles, wall_particles, smooth_kernel, wc_lm, dr_fti_forward_euler, vtp_io; check = check);
