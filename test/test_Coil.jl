using CoilForces
using Test

@testset "Test winding pack angle functions" begin
    @testset "FrameCircle should return normal, binormal" begin
        curve = CurveCircle(2.3)
        frame = FrameCircle()
        n = 10
        for j in 1:n
            ϕ = (j - 1) * 2π / n
            differential_arclength, curvature, torsion, position, tangent, normal, binormal = Frenet_frame(curve, ϕ)
            p, q = get_frame(frame, ϕ, position, tangent, normal)
            @test p ≈ normal
            @test q ≈ binormal
            @test dot(tangent, cross(p, q)) ≈ 1
        end
    end

    @testset "FrameCentroid should return the correct vectors for a circle" begin
        curve = CurveCircle(2.3)
        frame = FrameCentroid(curve)
        n = 10
        for j in 1:n
            ϕ = (j - 1) * 2π / n
            differential_arclength, curvature, torsion, position, tangent, normal, binormal = Frenet_frame(curve, ϕ)
            p, q = get_frame(frame, ϕ, position, tangent, normal)
            @test p ≈ -normal
            @test q ≈ -binormal
            @test dot(tangent, cross(p, q)) ≈ 1
        end
    end

    @testset "FrameCentroid should return vectors that, with the tangent, form an orthonormal basis" begin
        curve = get_curve("hsx", 1)
        frame = FrameCentroid(curve)
        n = 10
        for j in 1:n
            ϕ = (j - 1) * 2π / n
            differential_arclength, curvature, torsion, position, tangent, normal, binormal = Frenet_frame(curve, ϕ)
            p, q = get_frame(frame, ϕ, position, tangent, normal)
            @test norm(tangent) ≈ 1
            @test norm(p) ≈ 1
            @test norm(q) ≈ 1
            @test dot(tangent, p) ≈ 0 atol=1e-15
            @test dot(tangent, q) ≈ 0 atol=1e-15
            @test dot(p, q) ≈ 0 atol=1e-15
            @test dot(tangent, cross(p, q)) ≈ 1
        end
    end

end