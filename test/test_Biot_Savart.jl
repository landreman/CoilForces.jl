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

    @testset "Test B from HSX coils at specified points" begin
        # Compare to reference values from simsopt computed by
        # 20221224-01-HSX_BiotSavart_simsopt_julia_benchmark

        r_eval = [1.42, 0.1, 0.04]

        current = -1.5e5

        data = [[-0.072104818545038  0.271757521790311  0.16363853189801 ]
            [-0.084720529664466  0.173067825353653  0.098354521790391]
            [-0.040274051646095  0.070919120085303  0.035942181167016]
            [-0.018659296051856  0.030788105746346  0.013455935703873]
            [-0.010230770055778  0.014060109565271  0.005680799518643]
            [-0.00611556754025   0.006710428046814  0.002653092435928]]

        for jcoil in 1:6
            curve = get_curve("hsx", jcoil)
            coil = Coil(curve, current, 0.0)
            B_fixed = B_filament_fixed(coil, r_eval, 1600)
            B_adaptive = B_filament_adaptive(coil, r_eval)
            B_simsopt = data[jcoil, :]
            #println("point 1, adaptive:", B_adaptive)
            #println("point 1, fixed vs adaptive:", maximum(abs.(B_fixed ./ B_adaptive .- 1)))
            #println("point 1, simsopt vs adaptive:", maximum(abs.(B_simsopt ./ B_adaptive .- 1)))
            @test maximum(abs.(B_fixed ./ B_adaptive .- 1)) < 1e-9
            @test maximum(abs.(B_adaptive ./ B_simsopt .- 1)) < 1e-12
        end
    end
end
