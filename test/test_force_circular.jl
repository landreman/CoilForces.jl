using CoilForces
using Test

# Compare to reference values from 
# "20221019-02 Convergence of julia calculation of circular coil self-force.docx"
# and
# 20221016_04_Finite_thickness_circular_coil.jl

@testset "Test force for a circular coil" begin
    @testset "Compare analytic solution to a reference value" begin
        # Major radius of coil [meters]
        R0 = 2.3

        # Minor radius of coil [meters]
        a = 0.1

        # Total current [Amperes]
        I = 3.1e6

        curve = CurveCircle(R0)
        coil = Coil(curve, I, a)
        @test analytic_force_per_unit_length(coil) ≈ 1.8655666361140157e6
    end

    @testset "Test that regularized Biot-Savart for a filament matches the analytic result" begin
        # Major radius of coil [meters]
        R0 = 2.3

        # Minor radius of coil [meters]
        a = 0.1

        # Total current [Amperes]
        I = 3.1e6

        curve = CurveCircle(R0)
        coil = Coil(curve, I, a)

        regularization = a * a / sqrt(exp(1))
        r_eval = [R0, 0, 0]
        B_fixed = B_filament_fixed(coil, r_eval, 1000, regularization=regularization)
        B_adaptive = B_filament_adaptive(coil, r_eval, regularization=regularization)
        @test B_fixed ≈ B_adaptive
        @test abs(B_adaptive[1]) < 1e-13
        @test abs(B_adaptive[2]) < 1e-13
        force_from_quadrature = I * B_adaptive[3]
        @test force_from_quadrature ≈ analytic_force_per_unit_length(coil) rtol=1e-3
    end

    @testset "Test that force calculation for a finite thickness coil matches the analytic result" begin
        # Major radius of coil [meters]
        R0 = 2.3

        # Minor radius of coil [meters]
        a = 0.23

        # Total current [Amperes]
        I = 3.1e6

        curve = CurveCircle(R0)
        coil = Coil(curve, I, a)

        reltol = 1e-3
        abstol = 1e-10

        ϕ = 0
        @time force = force_finite_thickness(coil, ϕ, reltol=reltol, abstol=abstol)
        @show force
        @test abs(force[2]) < 1e-13
        @test abs(force[3]) < 1e-8
        @test force[1] ≈ analytic_force_per_unit_length(coil) rtol=3e-3
    end

    @testset "Check interpolated hifi calculation for several aspect ratios" begin
        # Major radius of coil [meters]
        R0 = 2.3

        # Total current [Amperes]
        I = 3.1e6

        curve = CurveCircle(R0)

        # Aspect ratio 1:
        coil = Coil(curve, I, R0)
        @test interpolated_force_per_unit_length(coil) / analytic_force_per_unit_length(coil) ≈ 0.8672096485513339
        # The reference value above is from the file
        # circular_coil_high_fidelity_over_analytic_force_rtol_1.0e-7_atol_1.0e-7_2023-02-17T02:57:13.639.dat

        # Aspect ratio 10:
        coil = Coil(curve, I, R0 / 10)
        @test interpolated_force_per_unit_length(coil) / analytic_force_per_unit_length(coil) ≈ 0.9987209678389167

        # Aspect ratio 100:
        coil = Coil(curve, I, R0 / 100)
        @test interpolated_force_per_unit_length(coil) / analytic_force_per_unit_length(coil) ≈ 0.9999873239392976

        # Aspect ratio 1e4:
        coil = Coil(curve, I, R0 / 1e4)
        @test interpolated_force_per_unit_length(coil) / analytic_force_per_unit_length(coil) ≈ 1


    end

end
