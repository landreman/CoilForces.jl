function plot_force_for_HSX()
    coil_num = 2
    curve = get_curve("hsx", coil_num)

    current = -1.5e5

    # minor radius of conductor:
    a = 0.02

    coil = Coil(curve, current, a)
    regularization = a * a / sqrt(ℯ)

    # Number of grid points on which to report the force.
    # (For computing the force, adaptive quadrature will be used.)
    nϕ = 200

    ϕ = (collect(1:nϕ) .- 1) * 2π / nϕ
    force_per_unit_length = zeros(nϕ, 3)
    for j in 1:nϕ
        B = B_filament_adaptive(coil, γ(curve, ϕ[j]), regularization=regularization)
        force_per_unit_length[j, :] = current * cross(tangent(curve, ϕ[j]), B)
    end

    plot(ϕ, force_per_unit_length[:, 1], label="x")
    plot!(ϕ, force_per_unit_length[:, 2], label="y")
    plot!(ϕ, force_per_unit_length[:, 3], label="x")
    xlabel!("Coil parameter ϕ")
    title!("Force per unit length on HSX coil $(coil_num) [N/m]")
end

function plot_integrand()
    coil_num = 2
    curve = get_curve("hsx", coil_num)
    # curve = CurveCircle(2.2)

    current = -1.5e5

    # minor radius of conductor:
    a = 0.001

    # point at which to evaluate the force:
    ϕ0 = 0.0

    coil = Coil(curve, current, a)
    δ = a * a / sqrt(ℯ)

    nϕ = 1000

    #ϕp = (collect(1:nϕ) .- 1) * 2π / nϕ .- π
    plot_width = 0.1
    ϕp = collect(range(ϕ0 - plot_width, ϕ0 + plot_width, length=nϕ))
    #ϕp = collect(range(-π, π, length=nϕ))
    integrand = zeros(nϕ, 3)
    smoothed_integrand = zeros(nϕ, 3)
    smoother_integrand = zeros(nϕ, 3)
    smoothest_integrand = zeros(nϕ, 3)
    r_eval = γ(curve, ϕ0)
    r_prime = dγdϕ(curve, ϕ0)
    r_prime_prime = d2γdϕ2(curve, ϕ0)
    r_prime_prime_prime = d3γdϕ3(curve, ϕ0)
    dot_factor = dot(r_prime, r_prime_prime)
    dℓdϕ_squared = normsq(r_prime)
    for j in 1:nϕ
        integrand[j, :] = d_B_d_ϕ(coil, ϕp[j], r_eval, regularization=δ)
        smoothed_integrand[j, :] = (integrand[j, :] + μ0 * current / (4π) * 0.5 * cross(r_prime_prime, r_prime)
            * (2 - 2 * cos(ϕp[j] - ϕ0)) / (((2 - 2 * cos(ϕp[j] - ϕ0)) * dℓdϕ_squared + δ) ^ 1.5))

        smoother_integrand[j, :] = (
            integrand[j, :] + μ0 * current / (4π) * 
            (0.5 * cross(r_prime_prime, r_prime) * (2 - 2 * cos(ϕp[j] - ϕ0))
            + (1.0/3) * cross(r_prime_prime_prime, r_prime) * (2 - 2 * cos(ϕp[j] - ϕ0)) * sin(ϕp[j] - ϕ0)
            ) 
            / (((2 - 2 * cos(ϕp[j] - ϕ0)) * dℓdϕ_squared + δ) ^ 1.5))

        #smoothed_integrand[j, :] = (integrand[j, :] + μ0 * current / (4π) * 0.5 * cross(r_prime_prime, r_prime)
        #    * (2 - 2 * cos(ϕp[j] - ϕ0)) / (((2 - 2 * cos(ϕp[j] - ϕ0)) * dℓdϕ_squared 
        #    + 0 * (2 - 2 * cos(ϕp[j] - ϕ0)) * sin(ϕp[j] - ϕ0) * dot_factor
        #    + δ) ^ 1.5))

        smoothest_integrand[j, :] = (integrand[j, :] + μ0 * current / (4π) * 
            (0.5 * cross(r_prime_prime, r_prime) * (2 - 2 * cos(ϕp[j] - ϕ0))
            + (1.0/3) * cross(r_prime_prime_prime, r_prime) * (2 - 2 * cos(ϕp[j] - ϕ0)) * sin(ϕp[j] - ϕ0)
            )
             / (((2 - 2 * cos(ϕp[j] - ϕ0)) * dℓdϕ_squared 
            + 1 * (2 - 2 * cos(ϕp[j] - ϕ0)) * sin(ϕp[j] - ϕ0) * dot_factor
            + δ) ^ 1.5))
    end

    plot()
    xyz = "xyz"
    #plot(ϕ, integrand[:, 1], label="x")
    #plot!(ϕ, integrand[:, 2], label="y")
    #plot!(ϕ, integrand[:, 3], label="z")
    for j in 1:3
        #plot!(ϕp, integrand[:, j], label=xyz[j:j])
        #plot!(ϕp, smoothed_integrand[:, j], label="smoothed, " * xyz[j:j])
        plot!(ϕp, smoother_integrand[:, j], label="smoother, " * xyz[j:j], ls=:dot)
        plot!(ϕp, smoothest_integrand[:, j], label="smoothest, " * xyz[j:j], ls=:dash)
    end
    xlabel!("Coil parameter ϕ")
    title!("Regularized Biot-Savart integrand for HSX coil $(coil_num) at ϕ=$(ϕ0)")

end

function plot_force_convergence_single()
    coil_num = 2

    # point at which to evaluate the force:
    ϕ0 = 0.0

    curve = get_curve("hsx", coil_num)
    # curve = CurveCircle(2.2)

    current = -1.5e5

    # minor radius of conductor:
    a = 0.001

    coil = Coil(curve, current, a)
    δ = a * a / sqrt(ℯ)

    # Generate numbers of quadrature points to try:
    nns = 40
    ns = [Int(round(10 ^ x)) for x in range(1.0, 4.0, length=nns)]

    r_eval = γ(curve, ϕ0)
    tangent0 = tangent(curve, ϕ0)
    force_per_unit_length = zeros(nns)
    force_per_unit_length_singularity_subtraction = zeros(nns)
    for jn in 1:nns
        B = B_filament_fixed(coil, r_eval, ns[jn], regularization=δ)
        force_per_unit_length[jn] = current * norm(cross(tangent0, B))

        B = B_singularity_subtraction_fixed(coil, ϕ0, ns[jn], δ)
        force_per_unit_length_singularity_subtraction[jn] = current * norm(cross(tangent0, B))
  end

    scatter(ns, force_per_unit_length, xscale=:log10, label="original")
    scatter!(ns, force_per_unit_length_singularity_subtraction, label="singularity subtraction")
    xlabel!("number of quadrature points")
    ylabel!("Force per unit length [N / m]")
end

function plot_force_convergence_grid()
    # points at which to evaluate the force:
    nϕ = 5
    ϕs = [(j - 1) * 2π / nϕ for j in 1:nϕ]

    current = -1.5e5

    # minor radius of conductor:
    a = 0.001
    δ = a * a / sqrt(ℯ)

    # Generate numbers of quadrature points to try:
    nns = 40
    ns = [Int(round(10 ^ x)) for x in range(1.0, 4.0, length=nns)]

    num_coils = 6
    layout = (length(ϕs), num_coils)
    plots = Array{Any, 2}(undef, num_coils, length(ϕs))
    style = (
        markershape = :circle,
        markersize = 2,
        markerstrokewidth = 0,
    )
    scalefontsizes()
    scalefontsizes(0.5)
    for coil_num in 1:num_coils
        curve = get_curve("hsx", coil_num)    
        coil = Coil(curve, current, a)
    
        for jϕ in 1:length(ϕs)
            ϕ0 = ϕs[jϕ]

            r_eval = γ(curve, ϕ0)
            tangent0 = tangent(curve, ϕ0)
            force_per_unit_length = zeros(nns)
            force_per_unit_length_singularity_subtraction = zeros(nns)
            for jn in 1:nns
                B = B_filament_fixed(coil, r_eval, ns[jn], regularization=δ)
                force_per_unit_length[jn] = current * norm(cross(tangent0, B))

                B = B_singularity_subtraction_fixed(coil, ϕ0, ns[jn])
                force_per_unit_length_singularity_subtraction[jn] = current * norm(cross(tangent0, B))
            end

            plots[coil_num, jϕ] = plot(ns, force_per_unit_length, xscale=:log10, label="original"; style...)
            plot!(ns, force_per_unit_length_singularity_subtraction, label="singularity subtraction"; style...)
            xlabel!("number of quadrature points")
            ylabel!("Force per unit length [N / m]")
            title = @sprintf "HSX coil %d, ϕ=%.2f" coil_num ϕ0
            title!(title)
            xticks!([1e1, 1e2, 1e3, 1e4])
        end
    end
    plot(plots..., layout=layout, dpi=100, size=(1100, 850))
    savefig("plot_force_convergence_grid.pdf")
end

function plot_B_near_thick_circular_coil()
    # Plot |B| inside and outside the coil.
    # This part takes tens of seconds to run.

    # Major radius of coil [meters]
    R0 = 2.3

    # Minor radius of coil [meters]
    aminor = 0.2

    # Total current [Amperes]
    I = 3.1e6

    reltol = 1e-3
    abstol = 1e-5
    nx = 20
    nz = 19

    curve = CurveCircle(R0)
    coil = Coil(curve, I, aminor)

    xplot = collect(range(R0 - 2 * aminor, R0 + 2 * aminor, length=nx))
    zplot = collect(range(- 1.5 * aminor, 1.5 * aminor, length=nz))
    modB = zeros(nz, nx)

    for jx in 1:nx
        for jz in 1:nz
            r_eval = [xplot[jx], 0, zplot[jz]]
            B = B_finite_thickness(coil, r_eval, reltol=reltol, abstol=abstol)
            modB[jz, jx] = sqrt(B[1]^2 + B[2]^2 + B[3]^2)
        end
    end

    contour(xplot, zplot, modB)
    title!("|B| [Tesla]")
    xlabel!("x [meters]")
    ylabel!("z [meters]")
    nθ = 150
    θplot = collect(range(0, 2π, length=nθ))
    plot!(R0 .+ aminor * cos.(θplot), aminor * sin.(θplot), linewidth=3, color=:black, label="coil edge")

end

function save_high_fidelity_force_for_circular_coil()
    reltol = 1e-5
    #abstol = reltol * 1e+0
    abstol = 1e-9

    #a_over_R = 10 .^ collect(((-4):(0.5):(-0.5)))
    a_over_R = 10 .^ collect(((-3):(0.5):(-3)))
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

        ϕ = 0
        time_data = @timed force = force_finite_thickness(coil, ϕ, reltol=reltol, abstol=abstol)
        analytic = analytic_force_per_unit_length(coil)
        high_fidelity_over_analytic_force[ja] = force[1] / analytic
        times[ja] = time_data.time
        println("  time: $(time_data.time)  analytic force: $(analytic)  (high fidelity force) / (analytic force): ", high_fidelity_over_analytic_force[ja])
    end
    @show high_fidelity_over_analytic_force

    """
    filename = "circular_coil_high_fidelity_over_analytic_force_rtol_$(reltol)_atol_$(abstol).dat"
    open(filename, "w") do file
        write(file, "a/R, (high fidelity force)/(analytic force), time\n")
        for ja in 1:length(a_over_R)
            write(file, "$(a_over_R[ja]), $(high_fidelity_over_analytic_force[ja]), $(times[ja])\n")
        end
    end
    """
end


function save_high_fidelity_force_for_circular_coil_vary_tol()
    #a_over_R = 10 .^ collect(((-4):(0.5):(-0.5)))
    a_over_R = 10 .^ collect(((-4):(0.5):(-4)))
    println("Values of a/R that will be evaluated: ", a_over_R)

    # Major radius of coil [meters]
    R0 = 1.0
    curve = CurveCircle(R0)

    # Total current [Amperes]
    I = 1.0

    high_fidelity_over_analytic_force = similar(a_over_R)
    times = similar(a_over_R)
    abstols = similar(a_over_R)
    reltols = similar(a_over_R)
    for ja in 1:length(a_over_R)
        a = a_over_R[ja]
        coil = Coil(curve, I, a)

        reltol = 1e-2
        #abstol = reltol * 1e+0
        abstol = 1e-12
        reltols[ja] = reltol
        abstols[js] = abstol
        println("a = $(a),  abstol = $(abstol),  reltol = $(reltol)")
        
        ϕ = 0
        time_data = @timed force = force_finite_thickness(coil, ϕ, reltol=reltol, abstol=abstol)
        analytic = analytic_force_per_unit_length(coil)
        high_fidelity_over_analytic_force[ja] = force[1] / analytic
        times[ja] = time_data.time
        println("  time: $(time_data.time)  analytic force: $(analytic)  (high fidelity force) / (analytic force): ", high_fidelity_over_analytic_force[ja])
    end
    @show high_fidelity_over_analytic_force

    """
    filename = "circular_coil_high_fidelity_over_analytic_force_rtol_$(reltol)_atol_$(abstol).dat"
    open(filename, "w") do file
        write(file, "a/R, (high fidelity force)/(analytic force), time\n")
        for ja in 1:length(a_over_R)
            write(file, "$(a_over_R[ja]), $(high_fidelity_over_analytic_force[ja]), $(times[ja])\n")
        end
    end
    """
end

function save_high_fidelity_force_for_circular_coil_tol_scan()
    a_over_R = 0.1

    # Major radius of coil [meters]
    R0 = 1.0
    curve = CurveCircle(R0)

    # Total current [Amperes]
    #I = 1.0
    I = 3.1e6

    coil = Coil(curve, I, a_over_R)

    reltols = 10 .^ collect(((-6):(1.0):(-1)))
    #reltols = [1e-2, 1e-4]
    println("Values of reltol that will be evaluated: ", reltols)
    n_reltols = length(reltols)

    abstols = 10 .^ collect(((-5):(1.0):(-1)))
    #abstols = [1e-6, 1e-4, 1e-2]
    println("Values of abstol that will be evaluated: ", abstols)
    n_abstols = length(abstols)

    high_fidelity_over_analytic_force = zeros(n_reltols, n_abstols)
    times = zeros(n_reltols, n_abstols)
    for ja in 1:n_abstols
        for jr in 1:n_reltols
            abstol = abstols[ja]
            reltol = reltols[jr]

            println("abstol = $(abstol),  reltol = $(reltol)")
            
            ϕ = 0
            time_data = @timed force = force_finite_thickness(coil, ϕ, reltol=reltol, abstol=abstol)
            analytic = analytic_force_per_unit_length(coil)
            high_fidelity_over_analytic_force[jr, ja] = force[1] / analytic
            times[jr, ja] = time_data.time
            println("  time: $(time_data.time)  analytic force: $(analytic)  (high fidelity force) / (analytic force): ", high_fidelity_over_analytic_force[ja])
        end
    end
    @show high_fidelity_over_analytic_force

    directory = "/Users/mattland/Box/work23/20230214-07_circular_coil_high_fidelity_over_analytic_convergence/"
    filename = "circular_coil_high_fidelity_over_analytic_force_a$(a_over_R)_$(Dates.now()).dat"
    open(directory * filename, "w") do file

        write(file, "abstols:\n")
        for ja in 1:n_abstols
            if ja > 1 write(file, ",") end
            write(file, "$(abstols[ja])")
        end
        write(file, "\n")

        write(file, "reltols:\n")
        for jr in 1:n_reltols
            if jr > 1 write(file, ",") end
            write(file, "$(reltols[jr])")
        end
        write(file, "\n")

        write(file, "(high fidelity force)/(analytic force):\n")
        for jr in 1:n_reltols
            for ja in 1:n_abstols
                if ja > 1 write(file, ",") end
                write(file, "$(high_fidelity_over_analytic_force[jr, ja])")
            end
            write(file, "\n")
        end

        write(file, "time:\n")
        for jr in 1:n_reltols
            for ja in 1:n_abstols
                if ja > 1 write(file, ",") end
                write(file, "$(times[jr, ja])")
            end
            write(file, "\n")
        end
    end
end