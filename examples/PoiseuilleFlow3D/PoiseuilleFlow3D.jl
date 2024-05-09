#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/05/09 15:20:52
  @ license: MIT
  @ description:
 =#

using EtherSPHCPU
using FLoops

const dim = 3;
const l = 1e-2;
const r0 = l / 2;
const dr = l / 40;
const gap = dr;
const influence_radius = 3.0 * dr;

smooth_kernel = SmoothKernel(influence_radius, dim, WendlandC2)

const nu = 1e-6;
const rho_0 = 1000.0;
const poiseuille_flow_theory_solutio_order = 10;
const f = 8e-5;
const ux_max = r0^2 / 4 * f / nu;

const x0 = 0.0;
const y0 = 0.0;
const z0 = 0.0;
# x as axis derection
# while y-o-z plate as the plane of the flow

function uxTheory(y::Float64, z::Float64)::Float64
    return f / 4 / nu * (r0^2 - y^2 - z^2)
end

const c_0 = 10 * ux_max;
const p_0 = 0.05 * rho_0 * c_0^2;
const mu_0 = nu * rho_0;
const gamma = 7;
const body_force_vec = [f, 0.0, 0.0];
wc_lm = CommonWeaklyCompressibleLiquidModel(rho_0, c_0, p_0, gamma, mu_0, body_force_vec);

const dt = 0.1 * influence_radius / c_0;
const total_time = 200.0
const total_step = total_time / dt |> round |> Int;
const output_interval = 10;
const density_reinitialized_interval = 5
const check_interval = 5;
dr_fti_forward_euler =
    DensityReinitializedFixedTimeIntervalForwardEuler(dt, total_step, output_interval, density_reinitialized_interval, check_interval);

const step_digit = div(total_step, output_interval) |> Int |> string |> length;
const file_name = "poiseuille_flow_3d";
const file_suffix = ".vtp";
const dir_path = "./examples/PoiseuilleFlow3D/PoiseuilleFlow3DData/";
const fileld_symbols = [:rho_, :p_, :c_];
const fileld_names = ["Density", "Pressure", "SoundSpeed"];
vtp_io = VTPIO(step_digit, file_name, file_suffix, dir_path, fileld_symbols, fileld_names);

function createCylinderRingParticles(particle_type, dx, x_0, y_0, z_0, inner_r, outer_r, cylinder_len; rho = NaN)
    particles = particle_type[]
    len_num = cylinder_len / dx |> round |> Int
    len_dx = cylinder_len / len_num
    r_num = (outer_r - inner_r) / dx |> round |> Int
    r_dx = (outer_r - inner_r) / r_num
    for i_len in 1:len_num
        x_here = x_0 + (i_len - 0.5) * len_dx
        for i_r in 1:r_num
            r_here = inner_r + (i_r - 0.5) * r_dx
            n_theta = 2 * pi * r_here / r_dx |> round |> Int
            if n_theta == 0
                n_theta = 1
            end
            d_theta = 2 * pi / n_theta
            dr_here = 2 * pi * r_here / n_theta
            for i_theta in 1:n_theta
                theta_here = i_theta * d_theta
                y_here = y_0 + r_here * cos(theta_here)
                z_here = z_0 + r_here * sin(theta_here)
                par = particle_type(dim)
                par.x_vec_ .= [x_here, y_here, z_here]
                if isnan(rho)
                else
                    par.rho_ = rho
                    par.mass_ = len_dx * dr_here * r_dx * par.rho_
                end
                par.gap_ = (len_dx + dr_here + r_dx) / dim
                push!(particles, par)
            end
        end
    end
    return particles
end

const pipe_len = 2 * l;

fluid_particles = createCylinderRingParticles(CommonLiquidParticle, dr, x0, y0, z0, 0.0, r0, pipe_len; rho = rho_0);
for index in eachindex(fluid_particles)
    fluid_particles[index].p_ = p_0
    fluid_particles[index].rho_ = rho_0
    fluid_particles[index].c_ = c_0
    fluid_particles[index].v_vec_[1] = uxTheory(fluid_particles[index].x_vec_[2], fluid_particles[index].x_vec_[3])
end

const wall_thick = 3
const wall_inner_r = r0
const wall_outer_r = r0 + wall_thick * dr
const wall_len = pipe_len + 4 * wall_thick * dr

wall_particles = createCylinderRingParticles(CompulsiveWallParticle, dr, x0 - 2 * wall_thick * dr, y0, z0, wall_inner_r, wall_outer_r, wall_len);
for index in eachindex(wall_particles)
    y = wall_particles[index].x_vec_[2]
    z = wall_particles[index].x_vec_[3]
    n_vec = [0.0, -y, -z] / sqrt(y^2 + z^2)
    wall_particles[index].normal_vec_ .= n_vec
end

let_thick = 3

inlet_particles = createCylinderRingParticles(CommonLiquidParticle, dr, x0 - let_thick * dr, y0, z0, 0.0, r0, let_thick * dr; rho = rho_0);
for index in eachindex(inlet_particles)
    inlet_particles[index].p_ = p_0
    inlet_particles[index].rho_ = rho_0
    inlet_particles[index].c_ = c_0
    inlet_particles[index].v_vec_[1] = uxTheory(inlet_particles[index].x_vec_[2], inlet_particles[index].x_vec_[3])
end

outlet_particles = createCylinderRingParticles(CommonLiquidParticle, dr, x0 + pipe_len, y0, z0, 0.0, r0, let_thick * dr; rho = rho_0);
for index in eachindex(outlet_particles)
    outlet_particles[index].p_ = p_0
    outlet_particles[index].rho_ = rho_0
    outlet_particles[index].c_ = c_0
    outlet_particles[index].v_vec_[1] = uxTheory(outlet_particles[index].x_vec_[2], outlet_particles[index].x_vec_[3])
end

inletIntoCalculationDomainCheck(p::AbstractParticle) = p.x_vec_[1] > x0
outletOutCalculationDomainCheck(p::AbstractParticle) = p.x_vec_[1] > x0 + pipe_len
outletOutCheck(p::AbstractParticle) = p.x_vec_[1] > x0 + pipe_len + let_thick * dr
inletIntoCalculationDomainModify(p::AbstractParticle) = begin
    p.x_vec_[1] -= let_thick * dr
end

function udf(; t::Float64 = 0.0)::Nothing
    inlet_task = @async begin
        @floop @simd for index in eachindex(inlet_particles)
            inlet_particles[index].v_vec_[1] = uxTheory(inlet_particles[index].x_vec_[2], inlet_particles[index].x_vec_[3])
        end
    end
    outlet_task = @async begin
        @floop @simd for index in eachindex(outlet_particles)
            outlet_particles[index].v_vec_[2] = 0.0
            outlet_particles[index].v_vec_[3] = 0.0
        end
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
