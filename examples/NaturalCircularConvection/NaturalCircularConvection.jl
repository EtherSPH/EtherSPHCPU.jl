#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/04/09 00:51:58
  @ license: MIT
  @ description:
 =#

using EtherSPHCPU

const prandtl_number = 0.706;
const rayleigh_number = 4.7e4;

# case 0 (standard case): ratio_io = 2.6, c_0 = 2.0, p_0 = 0.025*..., rayleigh_number = 4.7e4
# case 1: ratio_io = 3., c_0 = 2.0, p_0 = 0.025*..., rayleigh_number = 1e3/4/5
# case 2: ratio_io = 6., c_0 = 2.0, p_0 = 0.025*..., rayleigh_number = 1e4/4/5
# case 3: ratio_io = 10., c_0 = 2.0, p_0 = 0.025*..., rayleigh_number = 1e3/4/5
# case 4: ratio_io = 2/3/5/10, c_0 = 2.0, p_0 = 0.025*..., rayleigh_number = 4.7e4

const r_outer = 1.0;
const ratio_io = 2.6;
const r_inner = r_outer / ratio_io;
const reference_length = r_outer - r_inner;

const g = 1.0;
const g_vec = [0.0, -g];
const beta = 0.05;
const mu_0 = 0.001;
const kappa = 1.0;

const t_outer = 0.0;
const t_inner = 1.0;
const delta_t = t_inner - t_outer;
const t_0 = t_outer;

const cp = prandtl_number * kappa / mu_0;
const rho_0 = sqrt(rayleigh_number * mu_0 * kappa / g / beta / (reference_length^3) / delta_t / cp);
const c_0 = 2.0;
const p_0 = 0.025 * rho_0 * c_0^2;
const gamma = 7;
const alpha = mu_0 / rho_0 / prandtl_number;

const th_wc_lm = ThermalWeaklyCompressibleLiquidModel(rho_0, c_0, p_0, gamma, mu_0, g_vec, t_0, kappa, cp, beta);

const dim = 2;
const dr = 0.02;
const gap = dr;
const influence_radius = 3.0 * dr;

const smooth_kernel = SmoothKernel(influence_radius, dim, WendlandC2);

const dt = 0.1 * influence_radius / c_0;
const total_time = 200.0;
const total_step = total_time / dt |> round |> Int;
const output_interval = 100;
const density_reinitialized_interval = 5;
const check_interval = 100;
const dr_fti_forward_euler =
    DensityReinitializedFixedTimeIntervalForwardEuler(dt, total_step, output_interval, density_reinitialized_interval, check_interval)

const step_digit = div(total_step, output_interval) |> Int |> string |> length;
const file_name = "natural_circular_convection";
const file_sufix = ".vtp";
const dir_path = "./examples/NaturalCircularConvection/NaturalCircularConvectionData"
const field_symbols = [:rho_, :p_, :c_, :t_];
const field_names = ["Density", "Pressure", "Soundspeed", "Temperature"];
vtp_io = VTPIO(step_digit, file_name, file_sufix, dir_path, field_symbols, field_names)

const x0 = 0.0;
const y0 = 0.0;

function createCircularLiquidParticles(r_start::RealType where {RealType <: AbstractFloat}, r_end::RealType where {RealType <: AbstractFloat})
    n_layer = ceil(Int, (r_end - r_start) / dr)
    dr_here = (r_end - r_start) / n_layer
    particles = ThermalLiquidParticle[]
    for i_layer in 1:n_layer
        r_layer = r_start + (i_layer - 0.5) * dr_here
        n_particles_layer = ceil(Int, 2 * pi * r_layer / dr_here)
        dr_layer = 2 * pi * r_layer / n_particles_layer
        for i_particle in 1:n_particles_layer
            d_theta = 2 * pi * i_particle / n_particles_layer
            particle = ThermalLiquidParticle(2)
            particle.x_vec_ = [x0 + r_layer * cos(d_theta), y0 + r_layer * sin(d_theta)]
            particle.rho_ = rho_0
            particle.p_ = p_0
            particle.c_ = c_0
            particle.t_ = t_0
            particle.mass_ = rho_0 * dr_layer * dr_here
            particle.gap_ = sqrt(dr_layer * dr_here)
            particle.kappa_ = kappa
            particle.cp_ = cp
            push!(particles, particle)
        end
    end
    return particles
end

fluid_particles = createCircularLiquidParticles(r_inner, r_outer);

function createCircularThermostaticWallParticles(
    r_start::RealType where {RealType <: AbstractFloat},
    n_layer::IntType where {IntType <: Integer},
    normal_vec_towards::IntType where {IntType <: Integer}, # -1 inner towards, 1 outer towards
    t::RealType where {RealType <: AbstractFloat},
)
    particles = ThermostaticWallParticle[]
    dr_here = dr
    for i_layer in 1:n_layer
        r_layer = r_start + (i_layer - 0.5) * dr_here
        n_particles_layer = ceil(Int, 2 * pi * r_layer / dr_here)
        dr_layer = 2 * pi * r_layer / n_particles_layer
        for i_particle in 1:n_particles_layer
            d_theta = 2 * pi * i_particle / n_particles_layer
            particle = ThermostaticWallParticle(2)
            particle.x_vec_ = [x0 + r_layer * cos(d_theta), y0 + r_layer * sin(d_theta)]
            particle.normal_vec_ = [cos(d_theta), sin(d_theta)] * normal_vec_towards
            particle.t_ = t
            particle.gap_ = sqrt(dr_layer * dr_here)
            particle.kappa_ = kappa
            push!(particles, particle)
        end
    end
    return particles
end

const wall_thick_number = 3;

inner_wall_particles = createCircularThermostaticWallParticles(r_inner - wall_thick_number * dr, wall_thick_number, 1, t_inner);
outer_wall_particles = createCircularThermostaticWallParticles(r_outer, wall_thick_number, -1, t_outer);

thermostatic_wall_particles = vcat(inner_wall_particles, outer_wall_particles);

main() = @fastmath solve!(fluid_particles, thermostatic_wall_particles, smooth_kernel, th_wc_lm, dr_fti_forward_euler, vtp_io)
