using CoilForces
using Test

@testset "Test winding pack angle functions" begin
    @testset "FrameCircle should return normal, binormal" begin
        curve = CurveCircle(2.3)
        frame = FrameCircle()
        n = 10
        for j in 1:n
            ϕ = (j - 1) * 2π / n
            p, q = get_frame(frame, ϕ)
            differential_arclength, curvature, torsion, dℓdϕ, tangent, normal, binormal = Frenet_frame(curve, ϕ)
            @test p ≈ normal
            @test q ≈ binormal
            @test dot(tangent, cross(p, q)) ≈ 1
        end
    end

    @testset "FrameCentroid should return the right vectors for a circle" begin
    end

    @testset "FrameCentroid should return vectors that, with the tangent, form an orthonormal basis" begin
    end

end