using CoilForces
using Test

@testset "Tests for inductance calculations" begin
    @testset "For a circular coil, regularized filament method should match analytic result" begin
        R0 = 1.7
        aminor = 0.01
        current = 1.4e5
        curve = CurveCircle(R0)
        coil = Coil(curve, current, aminor)

        @test analytic_inductance_for_circular_coil(coil) ≈ inductance_filament_adaptive(coil) rtol=1e-5
    end
end