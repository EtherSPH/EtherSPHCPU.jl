#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/09 20:36:57
  @ license: MIT
  @ description:
 =#

@testset "Particle" begin
    @testset "FixedParticle" begin
        @test CompulsiveWallParticle(2).gap_ == 0.0
        @test StaticVelocityWallParticle(2).v_vec_ == [0.0, 0.0]
        @test ThermostaticWallParticle(2).t_ == 0.0
    end
    @testset "MovableParticle" begin
        @test CommonLiquidParticle(2).mass_ == 0.0
        @test ThermalLiquidParticle(2).t_ == 0.0
    end
end
