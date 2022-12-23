using CoilForces
using Test

@testset "Test CurveXYZFourier basics" begin
    @testset "Test CurveXYZFourier constructor" begin
        xc = [1.1, 0.2]
        xs = [0.0, 3.1]
        yc = [0.9, -0.1]
        ys = [0.0, -0.3]
        zc = [-0.2, 1.4]
        zs = [-1.3, 0.2]
        c = CurveXYZFourier(xc, xs, yc, ys, zc, zs)
        @test c.n == 2
    end
end

@testset "Test Frenet frame functions" begin
    @testset "Test CurveCircle" begin
        R0 = 3.7
        c = CurveCircle(R0)
        ϕ = -0.3
        dℓdϕ, κ, τ = curvature_torsion(c, ϕ)
        @test dℓdϕ ≈ R0
        @test κ ≈ 1 / R0
        @test τ ≈ 0.0
    end
end
