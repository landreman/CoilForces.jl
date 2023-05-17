using CoilForces
using Test

@testset "Vector potential tests" begin
    @testset "For a circular coil, compare finite thickness A to analytic result" begin
        R0 = 1.3
        aminor = 0.001
        current = 1.7e5
        curve = CurveCircle(R0)
        coil = Coil(curve, current, aminor)

        nρ = 5
        nθ = 7
        for jρ in 1:nρ
            for jθ in 1:nθ
                ρ = (jρ - 1) / (nρ - 1)
                θ = 2π * (jθ - 1) / (nθ - 1)
                r_eval = [R0 + aminor * ρ * cos(θ), 0, aminor * ρ * sin(θ)]
                A_hifi = CoilForces.A_finite_thickness(
                    coil, 
                    r_eval;
                    reltol=1e-5,
                    abstol=1e-12,
                )

                # For this analytic result, see the end of section 7 in
                # 20230517-01 Self-inductance of a general-shape coil with circular cross-section.pdf
                A_analytic = μ0 * current / (4π) * (-ρ^2 - 3 + 2 * log(8 * R0 / aminor))
                @test abs(A_hifi[1]) < 1e-8
                @test abs(A_hifi[3]) < 1e-8
                @test A_hifi[2] ≈ A_analytic rtol=5e-4
            end
        end
    end
end
