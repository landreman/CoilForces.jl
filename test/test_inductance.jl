using CoilForces
using Test

@testset "Tests for inductance calculations" begin
    @testset "For a circular coil, regularized filament method should match analytic result" begin
        R0 = 1.7
        aminor = 0.01
        current = 1.4e5
        curve = CurveCircle(R0)
        coil = CoilCircularXSection(curve, current, aminor)

        @test analytic_inductance_for_circular_coil(coil) ≈ inductance_filament_adaptive(coil) rtol=1e-5
    end

    @testset "For a circular coil with circular x-section, high fidelity method should match analytic result" begin
        R0 = 1.7
        aminor = 0.01
        current = 1.4e5
        curve = CurveCircle(R0)
        coil = CoilCircularXSection(curve, current, aminor)

        L_analytic = analytic_inductance_for_circular_coil(coil)
        tol = 1e-4
        L_hifi = inductance_finite_thickness(coil, reltol=tol, abstol=tol)
        @test L_analytic ≈ L_hifi rtol=1e-5
        
    end

    @testset "For HSX, hifi method should match regularized filament method" begin
        aminor = 0.01
        current = 1.4e5
        curve = get_curve("hsx", 1)
        coil = CoilCircularXSection(curve, current, aminor)

        L_filament = inductance_filament_adaptive(coil)
        tol = 1e-3
        L_hifi = inductance_finite_thickness(coil, reltol=tol, abstol=tol)
        @test L_filament ≈ L_hifi rtol=3e-4
    end

    @testset "Analytic inductance formula for rectancular x-section should be unchanged if we swap a and b" begin
        curve = CurveCircle(3.1)
        current = 1.2e5
        n_ratio = 10
        d = 0.01  # Geometric mean of a and b
        for j_ratio in 1:n_ratio
            log_ratio = ((j_ratio - 1) / (n_ratio - 1) - 0.5) * 2
            ratio = 10.0 ^ log_ratio
            a = d * ratio
            b = d / ratio
            coil1 = CoilRectangularXSection(curve, current, a, b, WindingPackAngleZero())
            coil2 = CoilRectangularXSection(curve, current, b, a, WindingPackAngleZero())
            @test analytic_inductance_for_circular_coil(coil1) ≈ analytic_inductance_for_circular_coil(coil2)
        end
    end

    @testset "Compare analytic inductance formula for rectancular x-section to Weinstein's form" begin
        R = 3.1
        curve = CurveCircle(R)
        current = 1.2e5
        n_ratio = 10
        d = 0.01  # Geometric mean of a and b
        for j_ratio in 1:n_ratio
            log_ratio = ((j_ratio - 1) / (n_ratio - 1) - 0.5) * 2
            ratio = 10.0 ^ log_ratio
            a = d * ratio
            b = d / ratio
            coil = CoilRectangularXSection(curve, current, a, b, WindingPackAngleZero())
            # Expression from Rosa
            x = b / a
            Weinstein_formula = (
                μ0 * R * (log(8 * R / a) 
                    + (1.0 / 12) 
                    - π * x / 3
                    + (- 0.5 + 1 / (12 * x * x)) * log(1 + x * x)
                    + x * x / 12 * log(1 + 1 / (x * x))
                    + (2.0 / 3) * (x - 1 / x) * atan(x)
            ))
            @test analytic_inductance_for_circular_coil(coil) ≈ Weinstein_formula
        end
    end

    @testset "For a circular coil with rectangular x-section, high fidelity method should match analytic result" begin
        R0 = 6.1
        curve = CurveCircle(R0)
        current = 1.2e6
        n_ratio = 5
        d = R0 * 0.01  # Geometric mean of a and b
        tol = 1e-3
        L_analytic = zeros(n_ratio)
        L_hifi = zeros(n_ratio)
        ratios = zeros(n_ratio)
        for j_ratio in 1:n_ratio
            log_ratio = ((j_ratio - 1) / (n_ratio - 1) - 0.5) * 1
            ratio = 10.0 ^ log_ratio
            ratios[j_ratio] = ratio
            @show ratio
            a = d * sqrt(ratio)
            b = d / sqrt(ratio)
            coil = CoilRectangularXSection(curve, current, a, b, WindingPackAngleZero())
            L_analytic[j_ratio] = analytic_inductance_for_circular_coil(coil)
            L_hifi[j_ratio] = inductance_finite_thickness(coil, reltol=tol, abstol=tol)
        end

        """
        plot(ratios, L_analytic, label="analytic")
        plot!(ratios, L_hifi, label="hi fi")
        #plot(ratios, L_analytic ./ L_hifi)
        xlabel!("a / b")
        title!("tol = $(tol)")
        """

        @show L_analytic
        @show L_hifi
        @test L_analytic ≈ L_hifi rtol=3e-3
    end
end