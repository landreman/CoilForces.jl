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
end
