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

@testset "Test special functions for rectangular cross-section" begin
    @testset "Check k for square cross-section" begin
        for a in [0.2, 1, 3.3]
            @test CoilForces.rectangular_xsection_k(a, a) ≈ 2.556493222766492
        end
    end

    @testset "Check δ for square cross-section" begin
        for a in [0.2, 1, 3.3]
            @test CoilForces.rectangular_xsection_δ(a, a) ≈ 0.19985294779417703
        end
    end

    @testset "k and δ for rectangular x-section should be unchanged if we swap a and b" begin
        d = 0.01  # Geometric mean of a and b
        for ratio in [0.1, 3.7]
            a = d * ratio
            b = d / ratio
            @test CoilForces.rectangular_xsection_k(a, b) ≈ CoilForces.rectangular_xsection_k(b, a)
            @test CoilForces.rectangular_xsection_δ(a, b) ≈ CoilForces.rectangular_xsection_δ(b, a)
        end
    end

    @testset "Check limits a >> b and b >> a" begin
        ratios = [1.1e6, 2.2e4, 3.5e5]
        xs = [0.2, 1.0, 7.3]
        for ratio in ratios
            for x in xs
                # a >> b
                b = x
                a = b * ratio
                @test CoilForces.rectangular_xsection_k(a, b) ≈ (7.0 / 6) + log(a / b) rtol=1e-3
                @test CoilForces.rectangular_xsection_δ(a, b) ≈ a / (b * exp(3)) rtol=1e-3

                # b >> a
                a = x
                b = ratio * a
                @test CoilForces.rectangular_xsection_k(a, b) ≈ (7.0 / 6) + log(b / a) rtol=1e-3
                @test CoilForces.rectangular_xsection_δ(a, b) ≈ b / (a * exp(3)) rtol=1e-3
            end
        end
    end
end