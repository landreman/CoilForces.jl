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

    @testset "For a circular coil, high fidelity method should match analytic result" begin
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
end