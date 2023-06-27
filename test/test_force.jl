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
        coil = CoilCircularXSection(curve, I, a)
        @test analytic_force_per_unit_length(coil) ≈ 1.8655666361140157e6
    end

    @testset "Test that regularized Biot-Savart for a filament matches the analytic result. Circular x-section" begin
        # Major radius of coil [meters]
        R0 = 2.3

        # Minor radius of coil [meters]
        a = 0.1

        # Total current [Amperes]
        I = 3.1e6

        curve = CurveCircle(R0)
        coil = CoilCircularXSection(curve, I, a)

        force_from_quadrature = force_filament_adaptive(coil, 0.0)
        @test abs(force_from_quadrature[2]) < 1e-13
        @test abs(force_from_quadrature[3]) < 1e-13

        regularization = a * a / sqrt(exp(1))
        r_eval = [R0, 0, 0]
        B_fixed = B_filament_fixed(coil, r_eval, 1000, regularization=regularization)
        B_adaptive = B_filament_adaptive(coil, r_eval, regularization=regularization)
        @test B_fixed ≈ B_adaptive
        @test abs(B_adaptive[1]) < 1e-13
        @test abs(B_adaptive[2]) < 1e-13
        force_from_quadrature2 = I * B_adaptive[3]
        @test force_from_quadrature[1] ≈ force_from_quadrature2
        @test force_from_quadrature[1] ≈ analytic_force_per_unit_length(coil) rtol=1e-3
    end

    @testset "Test that regularized Biot-Savart for a filament matches the analytic result. Rectangular x-section" begin
        # Major radius of coil [meters]
        R0 = 2.3

        # Total current [Amperes]
        I = 3.1e6

        as = [0.01, 0.02, 0.04]

        b_over_as = 10 .^ collect(-1:0.25:1)

        curve = CurveCircle(R0)
        for a in as
            for b_over_a in b_over_as
                b = a * b_over_a
                coil = CoilRectangularXSection(curve, I, a, b, FrameCircle())

                force_from_quadrature = force_filament_adaptive(coil, 0.0)
                @test abs(force_from_quadrature[2]) < 1e-13
                @test abs(force_from_quadrature[3]) < 1e-13

                regularization = CoilForces.compute_regularization(coil)
                r_eval = [R0, 0, 0]
                B_fixed = B_filament_fixed(coil, r_eval, 30000, regularization=regularization)
                B_adaptive = B_filament_adaptive(coil, r_eval, regularization=regularization)
                @test B_fixed ≈ B_adaptive
                @test abs(B_adaptive[1]) < 1e-13
                @test abs(B_adaptive[2]) < 1e-13
                force_from_quadrature2 = I * B_adaptive[3]
                @test force_from_quadrature[1] ≈ force_from_quadrature2
                force_analytic = analytic_force_per_unit_length(coil)
                @test force_from_quadrature[1] ≈ force_analytic rtol=1e-3
                #println("a: $(a)  b: $(b)")
                #println("  quad:     $(force_from_quadrature)")
                #println("  analytic: $(force_analytic)   diff: $(force_analytic - force_from_quadrature)")
            end
        end
    end

    @testset "Test that force calculation for a finite thickness coil matches the analytic result. Circular xsection" begin
        # Major radius of coil [meters]
        R0 = 2.3

        # Minor radius of coil [meters]
        a = 0.23

        # Total current [Amperes]
        I = 3.1e6

        curve = CurveCircle(R0)
        coil = CoilCircularXSection(curve, I, a)

        reltol = 1e-3
        abstol = 1e-10

        ϕ = 0
        @time force = force_finite_thickness(coil, ϕ, reltol=reltol, abstol=abstol)
        @show force
        @test abs(force[2]) < 1e-13
        @test abs(force[3]) < 1e-8
        @test force[1] ≈ analytic_force_per_unit_length(coil) rtol=3e-3
    end

    @testset "Test that force calculation for a finite thickness coil matches the analytic result. Rectangular xsection" begin
        # Major radius of coil [meters]
        R0 = 2.3

        # Dimensions of the cross-section [meters]
        a = 0.23
        b = 0.4

        # Total current [Amperes]
        I = 3.1e6

        curve = CurveCircle(R0)
        coil = CoilRectangularXSection(curve, I, a, b, FrameCircle())

        reltol = 3e-3
        abstol = 1e-10

        ϕ = 0
        @time force = force_finite_thickness(coil, ϕ, reltol=reltol, abstol=abstol)
        force_analytic = analytic_force_per_unit_length(coil)
        println("numerical force: $(force)")
        println("analytic force:   $(force_analytic)   rel diff: $((force[1] - force_analytic) / force_analytic)")
        @test abs(force[2]) < 1e-13
        @test abs(force[3]) < 1e-8
        @test force[1] ≈ force_analytic rtol=1e-3
    end

    @testset "Test that force calculation for a finite thickness coil (using single 5D integral) matches the analytic result" begin
        # Major radius of coil [meters]
        R0 = 2.3

        # Minor radius of coil [meters]
        a = 0.23

        # Total current [Amperes]
        I = 3.1e6

        curve = CurveCircle(R0)
        coil = CoilCircularXSection(curve, I, a)

        reltol = 5e-2
        abstol = 1e-10

        ϕ = 0
        @time force = force_finite_thickness_5D(coil, ϕ, reltol=reltol, abstol=abstol)
        @show force
        @test abs(force[2]) < 1e-13
        @test abs(force[3]) < 1e-0
        @test force[1] ≈ analytic_force_per_unit_length(coil) rtol=3e-3
    end

    @testset "Test that force calculation for a finite thickness coil (with Siena's best-fit circle trick) matches the analytic result" begin
        # Major radius of coil [meters]
        R0 = 2.3

        # Minor radius of coil [meters]
        a = 0.23

        # Total current [Amperes]
        I = 3.1e6

        curve = CurveCircle(R0)
        coil = CoilCircularXSection(curve, I, a)

        reltol = 1e-3
        abstol = 1e-10

        ϕ = 0
        @time integral, force_from_best_fit_circle, force = force_finite_thickness_singularity_subtraction(coil, ϕ, reltol=reltol, abstol=abstol)
        @test integral ≈ [0, 0, 0]
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
        coil = CoilCircularXSection(curve, I, R0)
        @test interpolated_force_per_unit_length(coil) / analytic_force_per_unit_length(coil) ≈ 0.8672096485513339
        # The reference value above is from the file
        # circular_coil_high_fidelity_over_analytic_force_rtol_1.0e-7_atol_1.0e-7_2023-02-17T02:57:13.639.dat

        # Aspect ratio 10:
        coil = CoilCircularXSection(curve, I, R0 / 10)
        @test interpolated_force_per_unit_length(coil) / analytic_force_per_unit_length(coil) ≈ 0.9987209678389167

        # Aspect ratio 100:
        coil = CoilCircularXSection(curve, I, R0 / 100)
        @test interpolated_force_per_unit_length(coil) / analytic_force_per_unit_length(coil) ≈ 0.9999873239392976

        # Aspect ratio 1e4:
        coil = CoilCircularXSection(curve, I, R0 / 1e4)
        @test interpolated_force_per_unit_length(coil) / analytic_force_per_unit_length(coil) ≈ 1


    end

end
