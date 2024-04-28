#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/04/26 16:27:43
  @ license: MIT
  @ description:
 =#

using EtherSPHCPU
using FLoops

const dim = 2;
const dr = 2e-4;
const gap = dr;
const influence_radius = 3.0 * dr;

smooth_kernel = SmoothKernel(influence_radius, dim, WendlandC2)

const l = 1e-2;
const nu = 1e-6;
const rho_0 = 1000.0;
const poiseuille_flow_theory_solutio_order = 10;
const f = 8e-5;
const ux_max = f / 2 / nu / 4 * l^2;

const x0 = 0.0;
const y0 = 0.0;

function uxTheory(y::Float64, t::Float64)::Float64
    y = y - y0
    ux = f / 2 / nu * y * (l - y)
    for n in 0:poiseuille_flow_theory_solutio_order
        ux -= 4 * f * l^2 / nu / pi^3 / (2 * n + 1)^3 * sin(pi * y / l * (2 * n + 1)) * exp(-(2 * n + 1)^2 * pi^2 * nu * t / l^2)
    end
    return ux
end

const c_0 = 10 * ux_max;
const p_0 = 0.02 * rho_0 * c_0^2;
const mu_0 = nu * rho_0;
const gamma = 7;
const body_force_vec = [f, 0.0];
wc_lm = CommonWeaklyCompressibleLiquidModel(rho_0, c_0, p_0, gamma, mu_0, body_force_vec);

const dt = 0.1 * influence_radius / c_0;
const total_time = 200.0
const total_step = total_time / dt |> round |> Int;
const output_interval = 100;
const density_reinitialized_interval = 5
const check_interval = 5;
dr_fti_forward_euler =
    DensityReinitializedFixedTimeIntervalForwardEuler(dt, total_step, output_interval, density_reinitialized_interval, check_interval);

const step_digit = div(total_step, output_interval) |> Int |> string |> length;
const file_name = "poiseuille_flow_2d";
const file_suffix = ".vtp";
const dir_path = "./examples/PoiseuilleFlow2D/PoiseuilleFlow2DData/";
const fileld_symbols = [:rho_, :p_, :c_];
const fileld_names = ["Density", "Pressure", "SoundSpeed"];
vtp_io = VTPIO(step_digit, file_name, file_suffix, dir_path, fileld_symbols, fileld_names);

function createRectangleParticles(partricle_type, dx, x_0, y_0, width_number, height_number)
    particles = [partricle_type(2) for _ in 1:(width_number * height_number)]
    for i_row in 1:height_number, j_col in 1:width_number
        index = (i_row - 1) * width_number + j_col
        particles[index].x_vec_ .= [x_0 + (j_col - 1) * dx, y_0 + (i_row - 1) * dx] .+ dx / 2
    end
    return particles
end

const fluid_width = l * 4;
const fluid_height = l;
const fluid_row_number = fluid_height / dr |> round |> Int;
const fluid_col_number = fluid_width / dr |> round |> Int;
const fluid_number = fluid_row_number * fluid_col_number;
fluid_particles = createRectangleParticles(CommonLiquidParticle, dr, x0, y0, fluid_col_number, fluid_row_number);
for index in 1:fluid_number
    fluid_particles[index].p_ = p_0
    fluid_particles[index].rho_ = rho_0
    fluid_particles[index].c_ = c_0
    fluid_particles[index].mass_ = dr^2 * rho_0
    fluid_particles[index].gap_ = gap
    fluid_particles[index].v_vec_[1] = uxTheory(fluid_particles[index].x_vec_[2], 0.0)
end

const inlet_thick_number = 3;
const inlet_width_number = inlet_thick_number;
const inlet_height_number = fluid_row_number;
const inlet_number = inlet_width_number * inlet_height_number;
inlet_particles = createRectangleParticles(CommonLiquidParticle, dr, x0 - inlet_thick_number * dr, y0, inlet_width_number, inlet_height_number);
for index in 1:inlet_number
    inlet_particles[index].p_ = p_0
    inlet_particles[index].rho_ = rho_0
    inlet_particles[index].c_ = c_0
    inlet_particles[index].mass_ = dr^2 * rho_0
    inlet_particles[index].gap_ = gap
    inlet_particles[index].v_vec_[1] = uxTheory(inlet_particles[index].x_vec_[2], 0.0)
end

const outlet_thick_number = 3;
const outlet_width_number = outlet_thick_number;
const outlet_height_number = fluid_row_number;
const outlet_number = outlet_width_number * outlet_height_number;
outlet_particles = createRectangleParticles(CommonLiquidParticle, dr, x0 + fluid_width, y0, outlet_width_number, outlet_height_number);
for index in 1:outlet_number
    outlet_particles[index].p_ = p_0
    outlet_particles[index].rho_ = rho_0
    outlet_particles[index].c_ = c_0
    outlet_particles[index].mass_ = dr^2 * rho_0
    outlet_particles[index].gap_ = gap
    outlet_particles[index].v_vec_[1] = uxTheory(outlet_particles[index].x_vec_[2], 0.0)
end

const wall_thick_number = 3;
const n_more = 3
const wall_wall_width_number = (fluid_width / dr |> round |> Int) + 2 * (wall_thick_number + n_more);
const wall_wall_height_number = wall_thick_number;

bottom_wall_particles = createRectangleParticles(
    CompulsiveWallParticle,
    dr,
    x0 - (wall_thick_number + n_more) * dr,
    y0 - wall_thick_number * dr,
    wall_wall_width_number,
    wall_thick_number,
);
for index in eachindex(bottom_wall_particles)
    bottom_wall_particles[index].normal_vec_ = [0.0, 1.0]
    bottom_wall_particles[index].gap_ = gap
end

top_wall_particles = createRectangleParticles(
    CompulsiveWallParticle,
    dr,
    x0 - (wall_thick_number + n_more) * dr,
    y0 + fluid_height,
    wall_wall_width_number,
    wall_thick_number,
);
for index in eachindex(top_wall_particles)
    top_wall_particles[index].normal_vec_ = [0.0, -1.0]
    top_wall_particles[index].gap_ = gap
end

wall_particles = vcat(bottom_wall_particles, top_wall_particles)

inletIntoCalculationDomainCheck(p::AbstractParticle) = p.x_vec_[1] > x0
outletOutCalculationDomainCheck(p::AbstractParticle) = p.x_vec_[1] > x0 + fluid_width
outletOutCheck(p::AbstractParticle) = p.x_vec_[1] > x0 + fluid_width + outlet_thick_number * dr
inletIntoCalculationDomainModify(p::AbstractParticle) = begin
    p.x_vec_[1] -= inlet_thick_number * dr
end

function udf(; t::Float64 = 0.0)::Nothing
    inlet_task = @async begin
        @floop @simd for index in eachindex(inlet_particles)
            inlet_particles[index].v_vec_[1] = uxTheory(inlet_particles[index].x_vec_[2], t)
        end
    end
    outlet_task = @async begin
        @floop @simd for index in eachindex(outlet_particles)
            outlet_particles[index].v_vec_[1] = uxTheory(outlet_particles[index].x_vec_[2], t)
            outlet_particles[index].v_vec_[2] = 0.0
        end
        nothing
    end
    wait(inlet_task)
    wait(outlet_task)
    return nothing
end

main() = @fastmath solve!(
    fluid_particles,
    inlet_particles,
    outlet_particles,
    wall_particles,
    smooth_kernel,
    wc_lm,
    dr_fti_forward_euler,
    vtp_io;
    inlet_into_calculation_domain_check = inletIntoCalculationDomainCheck,
    outlet_out_calculation_domain_check = outletOutCalculationDomainCheck,
    outlet_out_check = outletOutCheck,
    inlet_into_calculation_domain_modify = inletIntoCalculationDomainModify,
    user_defined_function = udf,
)
