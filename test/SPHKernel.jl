#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/09 17:31:46
  @ license: MIT
  @ description:
 =#

@testset "SPHKernel" begin
    @testset "SmoothKernel" begin
        dim = 2
        influence_radius = 1.0
        cubic_spline = SmoothKernel(influence_radius, dim, CubicSpline)
        @test isapprox(cubic_spline.h_, influence_radius / 2.0)
        gaussian = SmoothKernel(influence_radius, dim, Gaussian)
        @test isapprox(gaussian.h_, influence_radius / 3.0)
        wendland_c2 = SmoothKernel(influence_radius, dim, WendlandC2)
        @test isapprox(wendland_c2.h_, influence_radius / 2.0)
        wendland_c4 = SmoothKernel(influence_radius, dim, WendlandC4)
        @test isapprox(wendland_c4.h_, influence_radius / 2.0)
    end
end
