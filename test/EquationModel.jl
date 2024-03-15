#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/09 22:08:19
  @ license: MIT
  @ description:
 =#

@testset "EquationModel" begin
    @testset "LiquidModel" begin
        rho_0 = 1.0
        c_0 = 1.0
        p_0 = 1.0
        gamma = 7
        mu_0 = 0.1
        body_force_vec = [0.0, 0.0, -9.8]
        wc_lm = CommonWeaklyCompressibleLiquidModel(rho_0, c_0, p_0, gamma, mu_0, body_force_vec)
        @test isapprox(wc_lm.nu_0_, mu_0 / rho_0)
        t_0 = 1.0
        kappa_0 = 1.0
        cp = 1.0
        beta = 1.0
        twc_lm = ThermalWeaklyCompressibleLiquidModel(rho_0, c_0, p_0, gamma, mu_0, body_force_vec, t_0, kappa_0, cp, beta)
        @test isapprox(twc_lm.b_, c_0^2 * rho_0 / gamma)
    end
end
