using EtherSPHCPU

const prandtl_number = 0.706;
const rayleigh_number = 4.7e4;

const r_0 = 1.0;
const R = 2.6;
const r_i = r_0 / R;
const L = r_0 - r_i;

const rho_0 = 1.0;
const c_0 = 3.0;
const nu_0 = prandtl_number / sqrt(rayleigh_number);
const mu_0 = rho_0 * nu_0;
const p_0 = 0.02 * rho_0 * c_0^2;
const gamma = 7;
const alpha = nu_0 / prandtl_number;

const g = 1.0;
const g_vec = [0.0, -g];

const T_0 = 0.0;
const T_i = 1.0;
const delta_T = T_i - T_0;
const t_0 = (T_0 + T_i) / 2;
const kappa = 1.0;
const cp = prandtl_number * kappa / mu_0;
const beta = rayleigh_number * nu_0 * alpha / g / delta_T / (L^3);

const th_wc_lm = ThermalWeaklyCompressibleLiquidModel(rho_0, c_0, p_0, gamma, mu_0, g_vec, t_0, kappa, cp, beta);

const dim = 2;
const dr = 0.02;
const gap = dr;
const influence_radius = 3.0 * dr;

const smooth_kernel = SmoothKernel(influence_radius, dim, WendlandC4);

const dt = 0.1 * influence_radius / c_0;
const total_time = 30.0;
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

function creatCircularHeatLiquidParticles(r_start::RealType where {RealType <: AbstractFloat}, r_end::RealType where {RealType <: AbstractFloat})
    n_layer = ceil(Int, (r_end - r_start) / dr)
    particles = ThermalLiquidParticle[]
    for i_layer in 1:n_layer
        r_layer = r_start + (i_layer - 1) * dr + dr / 2
        n_particles_layer = ceil(Int, 2 * pi * r_layer / dr)
        dr_layer = 2 * pi * r_layer / n_particles_layer
        for i_particle in 1:n_particles_layer
            d_theta = 2 * pi * i_particle / n_particles_layer
            particle = ThermalLiquidParticle(2)
            particle.x_vec_ .= [x0 + r_layer * cos(d_theta), y0 + r_layer * sin(d_theta)]
            particle.rho_ = rho_0
            particle.p_ = p_0
            particle.c_ = c_0
            particle.t_ = t_0
            particle.mass_ = (dr_layer^2) * particle.rho_
            particle.gap_ = dr_layer
            particle.kappa_ = kappa
            push!(particles, particle)
        end
    end
    return particles
end

fluid_particles = creatCircularHeatLiquidParticles(r_i, r_0);

function creatCircularHeatSolidParticles(
    r_start::RealType where {RealType <: AbstractFloat},
    n_layer::IntType where {IntType <: Integer},
    normal_vec_towards::IntType where {IntType <: Integer}, # -1 inner towards, 1 outer towards
    t::RealType where {RealType <: AbstractFloat},
)
    particles = ThermostaticWallParticle[]
    for i_layer in 1:n_layer
        r_layer = r_start + (i_layer - 1) * dr + dr / 2
        n_particles_layer = ceil(Int, 2 * pi * r_layer / dr)
        dr_layer = 2 * pi * r_layer / n_particles_layer
        for i_particle in 1:n_particles_layer
            d_theta = 2 * pi * i_particle / n_particles_layer
            particle = ThermostaticWallParticle(2)
            particle.x_vec_ .= [x0 + r_layer * cos(d_theta), y0 + r_layer * sin(d_theta)]
            particle.normal_vec_ .= [cos(d_theta), sin(d_theta)] * normal_vec_towards
            particle.t_ = t
            particle.gap_ = dr_layer
            particle.kappa_ = kappa
            push!(particles, particle)
        end
    end
    return particles
end

const cavity_thick_number = 3;

In_Circular_particles = creatCircularHeatSolidParticles(r_i - 3 * dr, 3, 1, T_i);
Out_Circular_particles = creatCircularHeatSolidParticles(r_0, 3, -1, T_0);

thermostat_particles = vcat(In_Circular_particles, Out_Circular_particles);

main() = @fastmath solve!(fluid_particles, thermostat_particles, smooth_kernel, th_wc_lm, dr_fti_forward_euler, vtp_io)
