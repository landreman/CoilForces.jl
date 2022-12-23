using CoilForces
using Test

@testset "Test Biot-Savart law" begin
    @testset "Test z axis of circular coil" begin
        # Major radius of coil [meters]
        R0 = 2.3

        # Total current [Amperes]
        I = 3.1e6

        curve = CurveCircle(R0)
        coil = Coil(curve, I, 0.0)

        nz = 100
        z = collect(range(-5, 5, length=nz))
        Bz_analytic = @. 0.5 * μ0 * I * R0^2 / ((R0^2 + z^2) ^ 1.5)
        B = zeros(3, nz)
        nϕ = 100
        for j in 1:nz
            r_eval = [0, 0, z[j]]
            B[:, j] = B_filament(coil, r_eval, nϕ)
        end
        # Bx and By should be 0:
        @test maximum(abs.(B[1, :])) < 1e-12
        @test maximum(abs.(B[2, :])) < 1e-12
        @test maximum(abs.(B[3, :] ./ Bz_analytic .- 1)) < 1e-12
        
    end
end
