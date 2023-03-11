using CoilForces
using Test

@testset "Test Biot-Savart law for high fidelity circular coil" begin
    @testset "For B along the z axis for a circular coil, compare to analytic formula." begin
        # Major radius of coil [meters]
        R0 = 2.3
        
        # Total current [Amperes]
        I = 3.1e6
        
        # Coil minor radius
        a = 0.001
        
        nz = 10
        z = collect(range(-5, 5, length=nz))
        Bz_analytic = @. 0.5 * μ0 * I * R0^2 / ((R0^2 + z^2) ^ 1.5)
        Bz_numerical = zeros(nz)
        nϕ = 100
        for j in 1:nz
            r_eval = [0, 0, z[j]]
            # Evaluate Bz at (x, y, z) = (0, 0, z[j]):
            Bz_numerical[j] = hifi_circular_coil_compute_Bz(R0, a, I, 0, z[j]; reltol=1e-11, abstol=1e-13)
        end
        @test maximum(abs.(Bz_numerical ./ Bz_analytic .- 1)) < 1e-7

    end

    @testset "Compare specialized routine for circular coil to the hifi routine for general curve shapes" begin
        # Major radius of coil [meters]
        R0 = 2.3

        # Minor radius of coil [meters]
        aminor = 0.5

        # Total current [Amperes]
        I = 3.1e6

        reltol = 1e-4
        abstol = 1e-10
        nx = 5
        nz = 3

        curve = CurveCircle(R0)
        coil = Coil(curve, I, aminor)

        xplot = collect(range(R0 - aminor, R0 + aminor, length=nx))
        zplot = collect(range(- aminor, aminor, length=nz))
        Bz_specialized = zeros(nz, nx)
        Bz_general = zeros(nz, nx)

        for jx in 1:nx
            for jz in 1:nz
                r_eval = [xplot[jx], 0, zplot[jz]]
                B = B_finite_thickness(coil, r_eval, reltol=reltol, abstol=abstol)
                Bz_general[jz, jx] = B[3]
                Bz_specialized[jz, jx] = hifi_circular_coil_compute_Bz(R0, aminor, I, xplot[jx], zplot[jz]; reltol=reltol, abstol=abstol)
            end
        end
        """
        print("Bz_specialized:")
        display(Bz_specialized)
        print("Bz_general:")
        display(Bz_general)
        print("Difference:")
        display(Bz_specialized - Bz_general)
        """

        @test Bz_specialized ≈ Bz_general rtol=1e-5
    end
end

@testset "Test force for a high fidelity circular coil" begin

    @testset "Test that force calculation for a finite thickness coil matches the analytic result" begin
        # Major radius of coil [meters]
        R0 = 2.3

        # Minor radius of coil [meters]
        a = 0.1

        # Total current [Amperes]
        I = 3.1e6

        curve = CurveCircle(R0)
        coil = Coil(curve, I, a)

        reltol = 1e-4
        abstol = 1e-9

        @time force = hifi_circular_coil_force(R0, a, I; reltol=reltol, abstol=abstol)
        @show force
        @test force ≈ analytic_force_per_unit_length(coil) rtol=1e-3
    end

    @testset "Test that the specialized and general high fidelity methods give the same result" begin
        # Major radius of coil [meters]
        R0 = 1.7

        # Minor radius of coil [meters]
        a = 1.4

        # Total current [Amperes]
        I = 3.1e6

        curve = CurveCircle(R0)
        coil = Coil(curve, I, a)

        @time force_specialized = hifi_circular_coil_force(R0, a, I; reltol=1e-4, abstol=1e-8)
        @show force_specialized
        @time force_general = force_finite_thickness(coil, 0, reltol=1e-3, abstol=1e-8)
        @show force_general
        @test force_specialized ≈ force_general[1] rtol=1e-3
        @time integral, _, force_general_singularity_subtraction = force_finite_thickness_singularity_subtraction(coil, 0, reltol=1e-3, abstol=1e-8)
        @show force_general_singularity_subtraction
        @test integral ≈ [0, 0, 0]
        @test force_specialized ≈ force_general_singularity_subtraction[1] rtol=1e-3
    end

    @testset "Compare high fidelity force calculation to reference values for low aspect ratio" begin
        reltol = 1e-4
        abstol = 1e-4

        a_over_R = 10 .^ collect(((-0.5):(0.125):(0)))
        println("Values of a/R that will be evaluated: ", a_over_R)

        # Major radius of coil [meters]
        R0 = 1.0
        curve = CurveCircle(R0)

        # Total current [Amperes]
        I = 1.0

        high_fidelity_over_analytic_force = similar(a_over_R)
        times = similar(a_over_R)
        for ja in 1:length(a_over_R)
            a = a_over_R[ja]
            println("a = ", a)
            coil = Coil(curve, I, a)

            time_data = @timed force = hifi_circular_coil_force(R0, a, I; reltol=reltol, abstol=abstol)
            analytic = analytic_force_per_unit_length(coil)
            high_fidelity_over_analytic_force[ja] = force / analytic
            times[ja] = time_data.time
            println("  time: $(time_data.time)  analytic force: $(analytic)  (high fidelity force) / (analytic force): ", high_fidelity_over_analytic_force[ja])
        end
        @show high_fidelity_over_analytic_force

        # Reference values from
        # ~/Box/work23/20230216-01_circular_coil_high_fidelity_over_analytic_convergence/circular_coil_high_fidelity_over_analytic_force_rtol_0.0001_atol_0.0001_2023-02-16T08:15:20.764.dat
        reference_values = [
            0.987053088134848,
            0.9768546133681183,
            0.958577885497475,
            0.9257715919286639,
            0.8672151058186467,
        ]
        @test high_fidelity_over_analytic_force ≈ reference_values
    end
end
