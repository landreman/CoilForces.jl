using CoilForces
using Test

@testset "Vector potential tests" begin
    @testset "For a circular coil with circular x-section, compare finite thickness A to analytic result" begin
        R0 = 1.3
        aminor = 0.001
        current = 1.7e5
        curve = CurveCircle(R0)
        coil = CoilCircularXSection(curve, current, aminor)

        nρ = 5
        nθ = 7
        for jρ in 1:nρ
            for jθ in 1:nθ
                ρ = (jρ - 1) / (nρ - 1)
                θ = 2π * (jθ - 1) / (nθ - 1)
                r_eval = [R0 - aminor * ρ * cos(θ), 0, aminor * ρ * sin(θ)]
                A_hifi = CoilForces.A_finite_thickness(
                    coil, 
                    r_eval;
                    reltol=1e-6,
                    abstol=1e-12,
                )

                # For the leading-order analytic result, see the end of section 7 in
                # 20230517-01 Self-inductance of a general-shape coil with circular cross-section.pdf
                # For the high-order analytic result, see eq (143) in
                # 20230528-01_B_field_for_general_coil_with_circular_cross_section_coil_using_vector_potential.pdf
                A_analytic = μ0 * current / (4π) * (
                    -ρ^2 - 3 + 2 * log(8 * R0 / aminor)
                    + (-3 - 0.25 * ρ^2 + log(8 * R0 / aminor)) * (aminor / R0) * ρ * cos(θ)
                )
                @test abs(A_hifi[1]) < 1e-8
                @test abs(A_hifi[3]) < 1e-8
                #@test A_hifi[2] ≈ A_analytic rtol=5e-4  # Looser rtol for high-order terms only
                @test A_hifi[2] ≈ A_analytic rtol=4e-7
            end
        end
    end

    @testset "For a circular coil with rectangular x-section, compare finite thickness A to analytic result" begin
        R0 = 1.3
        a = 0.002
        b = 0.003
        current = 1.7e5
        curve = CurveCircle(R0)
        frame = FrameCircle()
        coil = CoilRectangularXSection(curve, current, a, b, frame)

        H(q, p) = 0.25 * (
            (a / b) * q * q * atan(b * p / (a * q))
            + (b / a) * p * p * atan(a * q / (b * p))
            + p * q * log(a * q * q / (4b) + b * p * p / (4a))
        )

        nu = 5
        nv = 7
        A_numerical = zeros(nu, nv)
        A_analytic = zeros(nu, nv)
        for ju in 1:nu
            for jv in 1:nv
                u = ((ju - 1) / (nu - 1) - 0.5) * 2 * 0.9999
                v = ((jv - 1) / (nv - 1) - 0.5) * 2 * 0.9999
                r_eval = [R0 - 0.5 * u * a, 0, 0.5 * v * b]
                A_hifi = CoilForces.A_finite_thickness(
                    coil, 
                    r_eval;
                    reltol=1e-6,
                    abstol=1e-12,
                )

                # For the leading-order analytic result, see
                # 20230531-01 Self-inductance for coils with rectangular cross-section.lyx
                temp = -1 + log(64 * R0 * R0 / (a * b))
                for su in [-1, 1]
                    for sv in [-1, 1]
                        temp -= su * sv * H(u + su, v + sv)
                    end
                end
                A_analytic[ju, jv] = μ0 * current / (4π) * temp
                A_numerical[ju, jv] = A_hifi[2]
                #println("$(A_hifi)")
                @test abs(A_hifi[1]) < 1e-8
                @test abs(A_hifi[3]) < 1e-8
                @test A_hifi[2] ≈ A_analytic[ju, jv] rtol=6e-4
            end
        end

        """
        scalefontsizes()
        scalefontsizes(0.5)

        p1 = contour(A_numerical, title="Numerical")
        p2 = contour(A_analytic, title="Analytic, leading")
        p3 = contour(A_numerical - A_analytic, title="num - leading")
        plot(p1, p2, p3,
            layout=(2, 2),
            size=(600, 600),
        )
        """
    end
end
