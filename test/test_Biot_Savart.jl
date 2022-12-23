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
            B_fixed = B_filament_fixed(coil, r_eval, nϕ)
            B_adaptive = B_filament_adaptive(coil, r_eval)
            #B_adaptive = B_filament_adaptive(coil, r_eval, reltol=1e-3, abstol=1e-5)
            # Fixed-grid quadrature and adaptive quadrature should agree:
            @test maximum(abs.(B_adaptive - B_fixed)) < 1e-12
            B[:, j] = B_adaptive
        end
        # Bx and By should be 0:
        @test maximum(abs.(B[1, :])) < 1e-12
        @test maximum(abs.(B[2, :])) < 1e-12
        @test maximum(abs.(B[3, :] ./ Bz_analytic .- 1)) < 1e-12
        
    end

    @testset "Test B from circular coil at specified points" begin
        # Compare a few points against the elliptic integral formula for an
        # infinitesmally thin coil, as implemented in simsopt.field.CircularCoil
        # 20221016_05_Benchmark_finite_thickness_circular_coil

        # Major radius of coil [meters]
        R0 = 2.3

        # Total current [Amperes]
        I = 3.1e6

        curve = CurveCircle(R0)
        coil = Coil(curve, I, 0.0)

        r_eval = [1.8, 0.7, 0.4]
        B_fixed = B_filament_fixed(coil, r_eval, 100)
        B_adaptive = B_filament_adaptive(coil, r_eval)
        B_simsopt = [0.797460697498886, 0.3101236045829,   1.210433050274526]
        #println("point 1, adaptive:", B_adaptive)
        #println("point 1, fixed vs adaptive:", maximum(abs.(B_fixed ./ B_adaptive .- 1)))
        #println("point 1, simsopt vs adaptive:", maximum(abs.(B_simsopt ./ B_adaptive .- 1)))
        @test maximum(abs.(B_fixed ./ B_adaptive .- 1)) < 1e-9
        @test maximum(abs.(B_adaptive ./ B_simsopt .- 1)) < 1e-13

        r_eval = [-3.5, -2.7, -1.4]
        B_fixed = B_filament_fixed(coil, r_eval, 100)
        B_adaptive = B_filament_adaptive(coil, r_eval)
        B_simsopt = [0.051493866798744,  0.039723840101888, -0.037540647636196]
        #println("point 2, adaptive:", B_adaptive)
        #println("point 2, fixed vs adaptive:", maximum(abs.(B_fixed ./ B_adaptive .- 1)))
        #println("point 2, simsopt vs adaptive:", maximum(abs.(B_simsopt ./ B_adaptive .- 1)))
        @test maximum(abs.(B_fixed ./ B_adaptive .- 1)) < 1e-12
        @test maximum(abs.(B_adaptive ./ B_simsopt .- 1)) < 1e-13
    end
end
