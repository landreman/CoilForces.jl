using CoilForces
using Test

@testset "Test Biot-Savart law" begin
    @testset "For B along the z axis for a circular coil, compare to analytic formula. Filament coil." begin
        # Major radius of coil [meters]
        R0 = 2.3

        # Total current [Amperes]
        I = 3.1e6

        curve = CurveCircle(R0)
        coils = [
            CoilCircularXSection(curve, I, 0.0)
            CoilRectangularXSection(curve, I, 0.0, 0.0, FrameCircle())
            CoilRectangularXSection(curve, I, 0.0, 0.0, FrameCentroid(curve))
        ]

        nz = 100
        z = collect(range(-5, 5, length=nz))
        Bz_analytic = @. 0.5 * μ0 * I * R0^2 / ((R0^2 + z^2) ^ 1.5)
        B = zeros(3, nz)
        nϕ = 100
        for coil in coils
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
        
    end

    @testset "For B along the z axis for a circular coil, compare to analytic formula. Finite thickness circular-cross-section coil." begin
        # Major radius of coil [meters]
        R0 = 2.3

        # Total current [Amperes]
        I = 3.1e6

        # Coil minor radius
        a = 0.001

        curve = CurveCircle(R0)
        coil = CoilCircularXSection(curve, I, a)

        nz = 10
        z = collect(range(-5, 5, length=nz))
        Bz_analytic = @. 0.5 * μ0 * I * R0^2 / ((R0^2 + z^2) ^ 1.5)
        B_numerical = zeros(3, nz)
        nϕ = 100
        for j in 1:nz
            r_eval = [0, 0, z[j]]
            B_numerical[:, j] = B_finite_thickness(coil, r_eval, reltol=1e-7, abstol=1e-8)
        end
        # Bx and By should be 0:
        @test maximum(abs.(B_numerical[1, :])) < 1e-10
        @test maximum(abs.(B_numerical[2, :])) < 1e-10
        @test maximum(abs.(B_numerical[3, :] ./ Bz_analytic .- 1)) < 1e-7
        
    end

    @testset "For B along the z axis for a circular coil, compare to analytic formula. Finite thickness rectangular-cross-section coil." begin
        # Major radius of coil [meters]
        R0 = 2.3

        # Total current [Amperes]
        I = 3.1e6

        # Coil cross-section dimensions
        a = 0.001
        b = 0.0017

        curve = CurveCircle(R0)
        coils = [
            CoilRectangularXSection(curve, I, a, b, FrameCircle()),
            CoilRectangularXSection(curve, I, a, b, FrameCentroid(curve)),
        ]

        nz = 10
        z = collect(range(-5, 5, length=nz))
        Bz_analytic = @. 0.5 * μ0 * I * R0^2 / ((R0^2 + z^2) ^ 1.5)
        B_numerical = zeros(3, nz)
        nϕ = 100
        for coil in coils
            for j in 1:nz
                r_eval = [0, 0, z[j]]
                B_numerical[:, j] = B_finite_thickness(coil, r_eval, reltol=1e-7, abstol=1e-8)
            end
            # Bx and By should be 0:
            @test maximum(abs.(B_numerical[1, :])) < 1e-10
            @test maximum(abs.(B_numerical[2, :])) < 1e-10
            @test maximum(abs.(B_numerical[3, :] ./ Bz_analytic .- 1)) < 1e-7
        end
        
    end

    @testset "For B from circular coil, compare to reference values from simsopt at specified points. Filament coil." begin
        # Compare a few points against the elliptic integral formula for an
        # infinitesmally thin coil, as implemented in simsopt.field.CircularCoil
        # 20221016_05_Benchmark_finite_thickness_circular_coil

        # Major radius of coil [meters]
        R0 = 2.3

        # Total current [Amperes]
        I = 3.1e6

        curve = CurveCircle(R0)
        coils = [
            CoilCircularXSection(curve, I, 0.0)
            CoilRectangularXSection(curve, I, 0.0, 0.0, FrameCircle())
            CoilRectangularXSection(curve, I, 0.0, 0.0, FrameCentroid(curve))
        ]

        for coil in coils
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

    @testset "For B from circular coil, compare to reference values from simsopt at specified points. Finite thickness circular-cross-section coil." begin
        # Compare a few points against the elliptic integral formula for an
        # infinitesmally thin coil, as implemented in simsopt.field.CircularCoil
        # 20221016_05_Benchmark_finite_thickness_circular_coil

        # Major radius of coil [meters]
        R0 = 2.3

        # Total current [Amperes]
        I = 3.1e6

        # Coil minor radius [meters]
        a = 0.001

        curve = CurveCircle(R0)
        coil = CoilCircularXSection(curve, I, a)

        r_eval = [1.8, 0.7, 0.4]
        B_julia = B_finite_thickness(coil, r_eval, reltol=1e-5, abstol=1e-7)
        B_simsopt = [0.797460697498886, 0.3101236045829,   1.210433050274526]
        #println("point 1, julia:  ", B_julia)
        #println("point 1, simsopt:", B_simsopt)
        #println("point 1, simsopt vs julia:", maximum(abs.(B_julia ./ B_simsopt .- 1)))
        @test maximum(abs.(B_julia ./ B_simsopt .- 1)) < 1e-6

        r_eval = [-3.5, -2.7, -1.4]
        B_julia = B_finite_thickness(coil, r_eval, reltol=1e-6, abstol=1e-7)
        B_simsopt = [0.051493866798744,  0.039723840101888, -0.037540647636196]
        #println("point 1, julia:  ", B_julia)
        #println("point 1, simsopt:", B_simsopt)
        #println("point 1, simsopt vs julia:", maximum(abs.(B_julia ./ B_simsopt .- 1)))
        @test maximum(abs.(B_julia ./ B_simsopt .- 1)) < 1e-7
    end

    @testset "For B from circular coil, compare to reference values from simsopt at specified points. Finite thickness rectangular-cross-section coil." begin
        # Compare a few points against the elliptic integral formula for an
        # infinitesmally thin coil, as implemented in simsopt.field.CircularCoil
        # 20221016_05_Benchmark_finite_thickness_circular_coil

        # Major radius of coil [meters]
        R0 = 2.3

        # Total current [Amperes]
        I = 3.1e6

        # Coil cross-section dimensions [meters]
        a = 0.0017
        b = 0.001

        curve = CurveCircle(R0)
        coils = [
            CoilRectangularXSection(curve, I, a, b, FrameCircle()),
            CoilRectangularXSection(curve, I, a, b, FrameCentroid(curve)),
        ]

        for coil in coils
            r_eval = [1.8, 0.7, 0.4]
            B_julia = B_finite_thickness(coil, r_eval, reltol=1e-5, abstol=1e-7)
            B_simsopt = [0.797460697498886, 0.3101236045829,   1.210433050274526]
            #println("point 1, julia:  ", B_julia)
            #println("point 1, simsopt:", B_simsopt)
            #println("point 1, simsopt vs julia:", maximum(abs.(B_julia ./ B_simsopt .- 1)))
            @test maximum(abs.(B_julia ./ B_simsopt .- 1)) < 1e-6

            r_eval = [-3.5, -2.7, -1.4]
            B_julia = B_finite_thickness(coil, r_eval, reltol=1e-6, abstol=1e-7)
            B_simsopt = [0.051493866798744,  0.039723840101888, -0.037540647636196]
            #println("point 1, julia:  ", B_julia)
            #println("point 1, simsopt:", B_simsopt)
            #println("point 1, simsopt vs julia:", maximum(abs.(B_julia ./ B_simsopt .- 1)))
            @test maximum(abs.(B_julia ./ B_simsopt .- 1)) < 2e-7
        end
    end

    @testset "For B from HSX coils, compare to reference values from simsopt at specified points" begin
        # Compare to reference values from simsopt computed by
        # 20221224-01-HSX_BiotSavart_simsopt_julia_benchmark

        r_eval = [1.42, 0.1, 0.04]

        current = -1.5e5

        # Coil minor radius, used only for the finite-thickness calculation
        aminor = 0.001

        data = [[-0.072104818545038  0.271757521790311  0.16363853189801 ]
            [-0.084720529664466  0.173067825353653  0.098354521790391]
            [-0.040274051646095  0.070919120085303  0.035942181167016]
            [-0.018659296051856  0.030788105746346  0.013455935703873]
            [-0.010230770055778  0.014060109565271  0.005680799518643]
            [-0.00611556754025   0.006710428046814  0.002653092435928]]

        for jcoil in 1:6
            curve = get_curve("hsx", jcoil)
            B_simsopt = data[jcoil, :]

            # First check circular x-section:
            coil = CoilCircularXSection(curve, current, aminor)
            B_fixed = B_filament_fixed(coil, r_eval, 1600)
            B_adaptive = B_filament_adaptive(coil, r_eval)
            B_thick = B_finite_thickness(coil, r_eval, reltol=1e-5, abstol=1e-7)
            #println("point 1, adaptive:", B_adaptive)
            #println("point 1, fixed vs adaptive:", maximum(abs.(B_fixed ./ B_adaptive .- 1)))
            #println("point 1, simsopt vs adaptive:", maximum(abs.(B_simsopt ./ B_adaptive .- 1)))
            @test maximum(abs.(B_fixed ./ B_adaptive .- 1)) < 1e-9
            @test maximum(abs.(B_adaptive ./ B_simsopt .- 1)) < 1e-12
            @test maximum(abs.(B_thick ./ B_simsopt .- 1)) < 1e-5

            # Now check rectangular x-section:
            coil = CoilRectangularXSection(curve, current, aminor, aminor * 1.3, FrameCentroid(curve))
            B_fixed = B_filament_fixed(coil, r_eval, 1600)
            B_adaptive = B_filament_adaptive(coil, r_eval)
            B_thick = B_finite_thickness(coil, r_eval, reltol=1e-5, abstol=1e-7)
            B_simsopt = data[jcoil, :]
            #println("point 1, adaptive:", B_adaptive)
            #println("point 1, fixed vs adaptive:", maximum(abs.(B_fixed ./ B_adaptive .- 1)))
            #println("point 1, simsopt vs adaptive:", maximum(abs.(B_simsopt ./ B_adaptive .- 1)))
            @test maximum(abs.(B_fixed ./ B_adaptive .- 1)) < 1e-9
            @test maximum(abs.(B_adaptive ./ B_simsopt .- 1)) < 1e-12
            @test maximum(abs.(B_thick ./ B_simsopt .- 1)) < 1e-5
        end
    end

    @testset "For a circular filament coil, check that the IxB force diverges logarithmically if the evaluation point is on the coil" begin
        # Major radius of coil [meters]
        R0 = 2.3

        # Total current [Amperes]
        I = 3.1e6

        curve = CurveCircle(R0)
        coil = CoilCircularXSection(curve, I, 0.0)
        r_eval = [R0, 0, 0]

        # Generate numbers of quadrature points to try:
        nns = 200
        ns = [Int(round(10 ^ x)) for x in range(2.0, 4.0, length=nns)]

        force_per_unit_length = zeros(nns)
        for jn in 1:nns
            B = B_filament_fixed(coil, r_eval, ns[jn], drop_first_point=true)
            force_per_unit_length[jn] = I * B[3]
        end

        # Fit a line in a semilog plot:
        c = (force_per_unit_length[end] - force_per_unit_length[1]) / (log(ns[end]) - log(ns[1]))
        d = force_per_unit_length[end] - c * log(ns[end])
        logarithmic_trend = c * log.(ns) .+ d   

        difference = force_per_unit_length ./ logarithmic_trend .- 1
        #@show difference
        @test maximum(abs.(difference)) < 2e-6

        if false
            using Plots
            scatter(ns, force_per_unit_length, xscale=:log10)
            plot!(ns, logarithmic_trend)
            xlabel!("number of quadrature points")
            ylabel!("Force per unit length [N / m]")
        end

    end

    @testset "Force from the singularity-subtraction method should match the force from direct quadrature of regularized Biot-Savart. Circular x-section." begin
        current = -1.5e5
    
        # minor radius of conductor:
        a = 0.001

        δ = a * a / sqrt(ℯ)

        # Number of points along each coil at which to evaluate the force:
        nϕ0 = 5

        # Number of points to use for quadrature:
        nϕ = 10000

        for coil_num in 1:6
            curve = get_curve("hsx", coil_num)
            coil = CoilCircularXSection(curve, current, a)

            # ϕ0 = point at which to evaluate the force:
            for ϕ0 in ((1:nϕ0) .- 1) * 2π / nϕ0
    
                r_eval, tangent0 = position_and_tangent(curve, ϕ0)

                B = B_filament_fixed(coil, r_eval, nϕ, regularization=δ)
                force_per_unit_length_original_fixed = current * cross(tangent0, B)
        
                B = B_filament_adaptive(coil, r_eval, regularization=δ)
                force_per_unit_length_original_adaptive = current * cross(tangent0, B)
        
                B = B_singularity_subtraction_fixed(coil, ϕ0, nϕ)
                force_per_unit_length_singularity_subtraction = current * cross(tangent0, B)

                @test force_per_unit_length_original_fixed ≈ force_per_unit_length_original_adaptive
                @test force_per_unit_length_original_adaptive ≈ force_per_unit_length_singularity_subtraction rtol=1e-5
            end
        end
    end

    @testset "Force from the singularity-subtraction method should match the force from direct quadrature of regularized Biot-Savart. Rectangular x-section." begin
        current = -1.5e5
    
        # cross-section dimensions:
        a = 0.001
        b = 0.0016

        regularization = a * b * CoilForces.rectangular_xsection_δ(a, b)

        # Number of points along each coil at which to evaluate the force:
        nϕ0 = 5

        # Number of points to use for quadrature:
        nϕ = 10000

        for coil_num in 1:6
            curve = get_curve("hsx", coil_num)
            coil = CoilRectangularXSection(curve, current, a, b, FrameCentroid(curve))

            # ϕ0 = point at which to evaluate the force:
            for ϕ0 in ((1:nϕ0) .- 1) * 2π / nϕ0
    
                r_eval, tangent0 = position_and_tangent(curve, ϕ0)

                B = B_filament_fixed(coil, r_eval, nϕ, regularization=regularization)
                force_per_unit_length_original_fixed = current * cross(tangent0, B)
        
                B = B_filament_adaptive(coil, r_eval, regularization=regularization)
                force_per_unit_length_original_adaptive = current * cross(tangent0, B)
        
                B = B_singularity_subtraction_fixed(coil, ϕ0, nϕ)
                force_per_unit_length_singularity_subtraction = current * cross(tangent0, B)

                @test force_per_unit_length_original_fixed ≈ force_per_unit_length_original_adaptive rtol=1e-6
                @test force_per_unit_length_original_adaptive ≈ force_per_unit_length_singularity_subtraction rtol=1e-6
            end
        end
    end

    @testset "Compare 2 versions of the finite-thickness Biot-Savart integrand" begin
        curve = get_curve("hsx", 2)
        coil = CoilCircularXSection(curve, 1.1e6, 0.05)
        ρ = 0.6
        θ = 0.3
        ϕ = 0.2
        r_eval = γ(curve, ϕ) + [0.001, 0.002, -0.003]
        @test CoilForces.B_finite_thickness_integrand(coil, ρ, θ, ϕ, r_eval) ≈ CoilForces.B_finite_thickness_integrand_sincos(coil, ρ, cos(θ), sin(θ), ϕ, r_eval)
    end

    @testset "For circular coil with circular x-section, compare hifi B vector to filament model" begin
        aminor = 0.01  # Minor radius of coil [meters]
        I = 1.0e4  # Total current [Amperes]
        reltol = 1e-4
        abstol = 1e-4
        nϕ = 1
        nρ = 4
        nθ = 3

        curve = CurveCircle(1.0)
        coil = CoilCircularXSection(curve, I, aminor)
        regularization = aminor * aminor / sqrt(ℯ)

        #####################################################
        #####################################################

        for jϕ in 1:nϕ
            #println("Processing jϕ = $(jϕ) of $(nϕ)")
            ϕ = 2π * (jϕ - 1) / nϕ
            differential_arclength, curvature, torsion, position, tangent, normal, binormal = Frenet_frame(curve, ϕ)
            B_regularized = B_filament_adaptive(coil, position; regularization=regularization)
            for jρ in 1:nρ
                #println("  Processing jρ = $(jρ) of $(nρ)")
                ρ = jρ / nρ
                for jθ = 1:nθ
                    θ = 2π * (jθ - 1) / nθ

                    # Evaluate vector B from a high fidelity calculation:
                    eval_point = position + aminor * ρ * cos(θ) * normal + aminor * ρ * sin(θ) * binormal
                    B_hifi = B_finite_thickness(coil, eval_point; reltol=1e-4, abstol=1e-4)

                    # Now evaluate vector B from the filament model:
                    B_filament = (
                        B_regularized
                        + CoilForces.B_local(coil, curvature, normal, binormal, ρ, θ)
                    )
                    #println("    B_hifi:     $(B_hifi)")
                    #println("    B_filament: $(B_filament)")
                    @test B_hifi ≈ B_filament atol=1e-12 rtol=3e-4
                end
            end
        end

    end

    @testset "For circular coil with rectangular x-section, compare hifi B vector to filament model" begin
        R0 = 2.0
        # Dimensions of the cross-section:
        a = 0.002
        b = 0.004
        I = 1.0e4  # Total current [Amperes]
        reltol = 1e-5
        abstol = 1e-5
        nx = 4
        nz = 3

        curve = CurveCircle(R0)
        coil = CoilRectangularXSection(curve, I, a, b, FrameCircle())
        r_eval = [R0, 0.0, 0.0]
        B_regularized_plus_extra_term = B_filament_adaptive(coil, r_eval; regularization=CoilForces.compute_regularization(coil))
        B_regularized_plus_extra_term[3] += μ0 * coil.current / (8π * R0) * (4 + 2 * log(2) + log(CoilForces.rectangular_xsection_δ(a, b)))
        κ1 = 1 / R0
        κ2 = 0.0

        #####################################################
        #####################################################

        almost_one = 0.999  # Avoid singularity at the coil edge
        for jx in 1:nx
            u = ((jx - 1) / (nx - 1) * 2 - 1) * almost_one
            x = R0 - 0.5 * u * a
            #println("jx = $(jx) of $(nx).  u = $(u)  x = $(x)")
            for jz = 1:nz
                v = ((jz - 1) / (nz - 1) * 2 - 1) * almost_one
                z = 0.5 * v * b
                #println("  jz = $(jz) of $(nz).  v = $(v)  z = $(z)")
    
                # Evaluate vector B from a high fidelity calculation:
                eval_point = [x, 0, z]
                B_hifi = B_finite_thickness(coil, eval_point; reltol=reltol, abstol=abstol)

                # Now evaluate vector B from the filament model:

                B0_p, B0_q = CoilForces.rectangular_xsection_B0(coil, u, v)
                B0 = [-B0_p, 0, B0_q]

                Bκ_p, Bκ_q = CoilForces.rectangular_xsection_Bκ(coil, κ1, κ2, u, v)
                Bκ_term = [-Bκ_p, 0, Bκ_q]

                B_filament = B_regularized_plus_extra_term + B0 + Bκ_term

                #println("    B_hifi:     $(B_hifi)")
                #println("    B_filament: $(B_filament)")
                #println("    diff:       $(B_hifi - B_filament)")
                @test B_hifi ≈ B_filament atol=1e-12 rtol=3e-6
            end
        end

    end
end
