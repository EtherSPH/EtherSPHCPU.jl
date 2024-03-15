#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/10 17:15:56
  @ license: MIT
  @ description:
 =#

@testset "Neighbour.jl" begin
    @test begin
        dim = 2
        influence_radius = 10.0
        smooth_kernel = SmoothKernel(influence_radius, dim, CubicSpline)
        par_1 = CommonLiquidParticle(dim)
        par_2 = CommonLiquidParticle(dim)
        par_1.x_vec_ .= [3.0, 4.0]
        neighbour_tuples = [(1, 1, 5.0)]
        neighbours = CommonNeighbour[]
        findNeighbours!(neighbour_tuples, neighbours, [par_1], [par_2], smooth_kernel)
        neighbour = neighbours[1]
        return isapprox(neighbour.r_, 5.0; atol = 1e-6)
    end
end
