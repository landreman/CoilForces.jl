"""
For a circular filament coil, show that the IxB force diverges logarithmically
if the evaluation point is on the coil and you merely skip the singular point.
"""
function plot_non_convergence_skipping_point_circular()
    # Major radius of coil [meters]
    R0 = 3.0

    # Minor radius of coil [meters]
    aminor = 0.01

    # Total current [Amperes]
    I = 1.0e6

    curve = CurveCircle(R0)
    coil = Coil(curve, I, aminor)
    r_eval = [R0, 0, 0]

    # Generate numbers of quadrature points to try:
    nns = 200
    ns = [Int(round(10 ^ x)) for x in range(1.0, 4.0, length=nns)]

    force_per_unit_length = zeros(nns)
    for jn in 1:nns
        B = B_filament_fixed(coil, r_eval, ns[jn], drop_first_point=true)
        force_per_unit_length[jn] = I * B[3]
    end

    x = [minimum(ns), maximum(ns)]

    linewidth = 2
    #scatter(ns, force_per_unit_length, xscale=:log10)
    plot(
        ns,
        force_per_unit_length,
        xscale=:log10,
        label="Skipping singular point",
        lw=linewidth,
        size=(600, 500),
        leg=false,
        c=:blue,
        xtickfont=font(10),
        ytickfont=font(10),
        xguidefontsize=12,
        yguidefontsize=12,
    )
    plot!(x, analytic_force_per_unit_length(coil) * [1, 1], label="Analytic", c=:red, lw=linewidth)
    plot!(x, interpolated_force_per_unit_length(coil) * [1, 1], linestyle=:dash, label="5D integral", c=:green, lw=4)
    xlabel!("Number of grid points for filament calculation")
    ylabel!("Force per unit length [N / m]")
    xticks!(10 .^ (1:4))
    xlims!(10, 1e4)
    ylims!((0, Inf))
    title!("Merely skipping the singular grid point for a filament         \ngives a non-convergent result with O(1) error        ")
    annotate!(25, 2.47e5, text("Analytic", :red))
    annotate!(30, 2.23e5, text("5D integral", :green))
    annotate!(180, 0.9e5, text("Filament, skipping singular point", :blue))
    savefig("circular_coil_force_non_convergence_skipping_point.pdf")
end

function plot_force_for_HSX()
    coil_num = 2
    curve = get_curve("hsx", coil_num)

    current = -1.5e5

    # minor radius of conductor:
    a = 0.02

    coil = Coil(curve, current, a)
    regularization = a * a / sqrt(???)

    # Number of grid points on which to report the force.
    # (For computing the force, adaptive quadrature will be used.)
    n?? = 200

    ?? = (collect(1:n??) .- 1) * 2?? / n??
    force_per_unit_length = zeros(n??, 3)
    for j in 1:n??
        B = B_filament_adaptive(coil, ??(curve, ??[j]), regularization=regularization)
        force_per_unit_length[j, :] = current * cross(tangent(curve, ??[j]), B)
    end

    plot(??, force_per_unit_length[:, 1], label="x")
    plot!(??, force_per_unit_length[:, 2], label="y")
    plot!(??, force_per_unit_length[:, 3], label="x")
    xlabel!("Coil parameter ??")
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
    ??0 = 0.0

    coil = Coil(curve, current, a)
    ?? = a * a / sqrt(???)

    n?? = 1000

    #??p = (collect(1:n??) .- 1) * 2?? / n?? .- ??
    plot_width = 0.1
    ??p = collect(range(??0 - plot_width, ??0 + plot_width, length=n??))
    #??p = collect(range(-??, ??, length=n??))
    integrand = zeros(n??, 3)
    smoothed_integrand = zeros(n??, 3)
    smoother_integrand = zeros(n??, 3)
    smoothest_integrand = zeros(n??, 3)
    r_eval = ??(curve, ??0)
    r_prime = d??d??(curve, ??0)
    r_prime_prime = d2??d??2(curve, ??0)
    r_prime_prime_prime = d3??d??3(curve, ??0)
    dot_factor = dot(r_prime, r_prime_prime)
    d???d??_squared = normsq(r_prime)
    for j in 1:n??
        integrand[j, :] = d_B_d_??(coil, ??p[j], r_eval, regularization=??)
        smoothed_integrand[j, :] = (integrand[j, :] + ??0 * current / (4??) * 0.5 * cross(r_prime_prime, r_prime)
            * (2 - 2 * cos(??p[j] - ??0)) / (((2 - 2 * cos(??p[j] - ??0)) * d???d??_squared + ??) ^ 1.5))

        smoother_integrand[j, :] = (
            integrand[j, :] + ??0 * current / (4??) * 
            (0.5 * cross(r_prime_prime, r_prime) * (2 - 2 * cos(??p[j] - ??0))
            + (1.0/3) * cross(r_prime_prime_prime, r_prime) * (2 - 2 * cos(??p[j] - ??0)) * sin(??p[j] - ??0)
            ) 
            / (((2 - 2 * cos(??p[j] - ??0)) * d???d??_squared + ??) ^ 1.5))

        #smoothed_integrand[j, :] = (integrand[j, :] + ??0 * current / (4??) * 0.5 * cross(r_prime_prime, r_prime)
        #    * (2 - 2 * cos(??p[j] - ??0)) / (((2 - 2 * cos(??p[j] - ??0)) * d???d??_squared 
        #    + 0 * (2 - 2 * cos(??p[j] - ??0)) * sin(??p[j] - ??0) * dot_factor
        #    + ??) ^ 1.5))

        smoothest_integrand[j, :] = (integrand[j, :] + ??0 * current / (4??) * 
            (0.5 * cross(r_prime_prime, r_prime) * (2 - 2 * cos(??p[j] - ??0))
            + (1.0/3) * cross(r_prime_prime_prime, r_prime) * (2 - 2 * cos(??p[j] - ??0)) * sin(??p[j] - ??0)
            )
             / (((2 - 2 * cos(??p[j] - ??0)) * d???d??_squared 
            + 1 * (2 - 2 * cos(??p[j] - ??0)) * sin(??p[j] - ??0) * dot_factor
            + ??) ^ 1.5))
    end

    plot()
    xyz = "xyz"
    #plot(??, integrand[:, 1], label="x")
    #plot!(??, integrand[:, 2], label="y")
    #plot!(??, integrand[:, 3], label="z")
    for j in 1:3
        #plot!(??p, integrand[:, j], label=xyz[j:j])
        #plot!(??p, smoothed_integrand[:, j], label="smoothed, " * xyz[j:j])
        plot!(??p, smoother_integrand[:, j], label="smoother, " * xyz[j:j], ls=:dot)
        plot!(??p, smoothest_integrand[:, j], label="smoothest, " * xyz[j:j], ls=:dash)
    end
    xlabel!("Coil parameter ??")
    title!("Regularized Biot-Savart integrand for HSX coil $(coil_num) at ??=$(??0)")

end

function plot_force_convergence_single()
    coil_num = 2

    # point at which to evaluate the force:
    ??0 = 0.0

    curve = get_curve("hsx", coil_num)
    # curve = CurveCircle(2.2)

    current = -1.5e5

    # minor radius of conductor:
    a = 0.001

    coil = Coil(curve, current, a)
    ?? = a * a / sqrt(???)

    # Generate numbers of quadrature points to try:
    nns = 40
    ns = [Int(round(10 ^ x)) for x in range(1.0, 4.0, length=nns)]

    r_eval = ??(curve, ??0)
    tangent0 = tangent(curve, ??0)
    force_per_unit_length = zeros(nns)
    force_per_unit_length_singularity_subtraction = zeros(nns)
    for jn in 1:nns
        B = B_filament_fixed(coil, r_eval, ns[jn], regularization=??)
        force_per_unit_length[jn] = current * norm(cross(tangent0, B))

        B = B_singularity_subtraction_fixed(coil, ??0, ns[jn], ??)
        force_per_unit_length_singularity_subtraction[jn] = current * norm(cross(tangent0, B))
  end

    scatter(ns, force_per_unit_length, xscale=:log10, label="original")
    scatter!(ns, force_per_unit_length_singularity_subtraction, label="singularity subtraction")
    xlabel!("number of quadrature points")
    ylabel!("Force per unit length [N / m]")
end

function plot_force_convergence_grid()
    # points at which to evaluate the force:
    n?? = 5
    ??s = [(j - 1) * 2?? / n?? for j in 1:n??]

    current = -1.5e5

    # minor radius of conductor:
    a = 0.001
    ?? = a * a / sqrt(???)

    # Generate numbers of quadrature points to try:
    nns = 40
    ns = [Int(round(10 ^ x)) for x in range(1.0, 4.0, length=nns)]

    num_coils = 6
    layout = (length(??s), num_coils)
    plots = Array{Any, 2}(undef, num_coils, length(??s))
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
    
        for j?? in 1:length(??s)
            ??0 = ??s[j??]

            r_eval = ??(curve, ??0)
            tangent0 = tangent(curve, ??0)
            force_per_unit_length = zeros(nns)
            force_per_unit_length_singularity_subtraction = zeros(nns)
            for jn in 1:nns
                B = B_filament_fixed(coil, r_eval, ns[jn], regularization=??)
                force_per_unit_length[jn] = current * norm(cross(tangent0, B))

                B = B_singularity_subtraction_fixed(coil, ??0, ns[jn])
                force_per_unit_length_singularity_subtraction[jn] = current * norm(cross(tangent0, B))
            end

            plots[coil_num, j??] = plot(ns, force_per_unit_length, xscale=:log10, label="original"; style...)
            plot!(ns, force_per_unit_length_singularity_subtraction, label="singularity subtraction"; style...)
            xlabel!("number of quadrature points")
            ylabel!("Force per unit length [N / m]")
            title = @sprintf "HSX coil %d, ??=%.2f" coil_num ??0
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
    n?? = 150
    ??plot = collect(range(0, 2??, length=n??))
    plot!(R0 .+ aminor * cos.(??plot), aminor * sin.(??plot), linewidth=3, color=:black, label="coil edge")

end

function save_high_fidelity_force_for_circular_coil_a_scan()
    #reltol = 3e-7
    #abstol = 3e-7
    reltol = 1e-6
    abstol = 1e-6

    #a_over_R = 10 .^ collect(((-4):(0.5):(-0.5)))
    #a_over_R = 10 .^ collect(((-2.25):(0.25):(-0.5)))
    #a_over_R = 10 .^ collect(((-1.0):(0.1):(-0.1)))
    a_over_R = 10 .^ collect(((-2.5):(0.0625):(0)))
    #a_over_R = 10 .^ collect(((0.0):(0.0625):(0.0)))
    println("Values of a/R that will be evaluated: ", a_over_R)

    # Major radius of coil [meters]
    R0 = 1.0
    curve = CurveCircle(R0)

    # Total current [Amperes]
    I = 1.0

    high_fidelity_over_analytic_force = similar(a_over_R)
    times = similar(a_over_R)
    for ja in 1:length(a_over_R)
        a = a_over_R[ja] * R0
        println("a = ", a)
        coil = Coil(curve, I, a)

        time_data = @timed force = hifi_circular_coil_force(R0, a, I; reltol=reltol, abstol=abstol)
        analytic = analytic_force_per_unit_length(coil)
        high_fidelity_over_analytic_force[ja] = force / analytic
        times[ja] = time_data.time
        println("  time: $(time_data.time)  analytic force: $(analytic)  (high fidelity force) / (analytic force): ", high_fidelity_over_analytic_force[ja])
    end
    @show high_fidelity_over_analytic_force

    #directory = "/Users/mattland/Box/work23/20230214-07_circular_coil_high_fidelity_over_analytic_convergence/"
    directory = "/Users/mattland/Box/work23/20230216-01_circular_coil_high_fidelity_over_analytic_convergence/"
    filename = "circular_coil_high_fidelity_over_analytic_force_rtol_$(reltol)_atol_$(abstol)_$(Dates.now()).dat"
    open(directory * filename, "w") do file
        write(file, "a/R, (high fidelity force)/(analytic force), time\n")
        for ja in 1:length(a_over_R)
            write(file, "$(a_over_R[ja]), $(high_fidelity_over_analytic_force[ja]), $(times[ja])\n")
        end
    end
    
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
        
        ?? = 0
        time_data = @timed force = force_finite_thickness(coil, ??, reltol=reltol, abstol=abstol)
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
            
            ?? = 0
            time_data = @timed force = force_finite_thickness(coil, ??, reltol=reltol, abstol=abstol)
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

function save_high_fidelity_force_for_circular_coil2_tol_scan()
    a_over_R = 0.01

    # Major radius of coil [meters]
    R0 = 1.0
    curve = CurveCircle(R0)

    # Total current [Amperes]
    I = 1.0
    
    coil = Coil(curve, I, a_over_R)

    reltols = 10 .^ collect(((-6):(1.0):(-2)))
    #reltols = [1e-2, 1e-4]
    println("Values of reltol that will be evaluated: ", reltols)
    n_reltols = length(reltols)

    abstols = 10 .^ collect(((-8):(1.0):(-3)))
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
            
            time_data = @timed force = hifi_circular_coil_force(R0, a_over_R, I; reltol=reltol, abstol=abstol)
            analytic = analytic_force_per_unit_length(coil)
            high_fidelity_over_analytic_force[jr, ja] = force / analytic
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

"""
This function reproduces a calculation from Siena's candidacy paper.
"""
function plot_Fx_for_HSX_coil_1_Siena()
    # Number of evaluation points for the force:
    neval = 300

    curve = get_curve("hsx", 1)
    current = 1.0
    aminor = 0.001
    coil = Coil(curve, current, aminor)
    regularization = aminor * aminor / sqrt(exp(1))

    force = zeros((neval, 3))
    force_circular_approximation = zeros((neval, 3))
    ?? = [2?? * (j - 1) / neval for j in 1:neval]
    for j in 1:neval
        position = ??(curve, ??[j])
        tangent_vector = tangent(curve, ??[j])
        B = B_filament_adaptive(coil, position; regularization=regularization)
        force[j, :] = current * cross(tangent_vector, B)
        force_circular_approximation[j, :] = force_locally_circular_approximation(coil, ??[j])
    end

    plot(??, force[:, 1], label="filament, adaptive quadrature")
    plot!(??, force_circular_approximation[:, 1], label="locally circular approximation")
    #plot!(??, force[:, 2], label="y")
    #plot!(??, force[:, 3], label="z")
    xlabel!("??")
    title!("dF???/d??? for HSX coil 1")
end

function plot_Fx_for_HSX_coil_1()
    # Number of evaluation points for the force:
    neval = 300

    curve = get_curve("hsx", 1)
    current = 1.0
    # Data from HsxCoilsNescToFinite.pdf:
    x_sectional_area = (56.8e-3) * (129.6e-3)
    aminor = sqrt(x_sectional_area / ??)
    @show aminor

    coil = Coil(curve, current, aminor)
    regularization = aminor * aminor / sqrt(exp(1))

    n_quadpoints_arr = [100, 1000, 10000]
    n_quadpointss = length(n_quadpoints_arr)
    force = zeros((neval, 3))
    force_circular_approximation = zeros((neval, 3))
    force_drop_singular_point = zeros((neval, 3, n_quadpointss))
    ?? = [2?? * (j - 1) / neval for j in 1:neval]
    for j in 1:neval
        position = ??(curve, ??[j])
        tangent_vector = tangent(curve, ??[j])
        B = B_filament_adaptive(coil, position; regularization=regularization)
        force[j, :] = current * cross(tangent_vector, B)
        force_circular_approximation[j, :] = force_locally_circular_approximation(coil, ??[j])
        for jquadpoints in 1:n_quadpointss
            n_quadpoints = n_quadpoints_arr[jquadpoints]
            d?? = 2?? / n_quadpoints
            B = [0.0, 0.0, 0.0]
            for k in 2:n_quadpoints
                ??_shifted = ??[j] + (k - 1) * d??
                #B += d_B_d_??(coil, ??_shifted, position, regularization=regularization)
                B += d_B_d_??(coil, ??_shifted, position, regularization=0)
            end
            B *= d??
            force_drop_singular_point[j, :, jquadpoints] = current * cross(tangent_vector, B)
        end
    end

    plot(??, force[:, 1], label="filament, adaptive quadrature")
    #plot!(??, force_circular_approximation[:, 1], label="locally circular approximation")
    for jquadpoints in 1:n_quadpointss
        plot!(??, force_drop_singular_point[:, 1, jquadpoints], label="filament, drop singular point, $(n_quadpoints_arr[jquadpoints]) points")
    end
    xlabel!("??")
    title!("dF???/d??? for HSX coil 1")
end


function save_high_fidelity_force_for_HSX()
    reltol = 1e-1
    abstol = 1e-1

    n?? = 4
    ??s = [2?? * (j - 1) / n?? for j in 1:n??]
    println("Values of ?? that will be evaluated: ", ??s)

    curve = get_curve("hsx", 1)
    #curve = CurveCircle(0.5)

    # Total current [Amperes]
    current = 150.0e3
    # Data from HsxCoilsNescToFinite.pdf:
    x_sectional_area = (56.8e-3) * (129.6e-3)
    aminor = sqrt(x_sectional_area / ??)
    coil = Coil(curve, current, aminor)


    forces = zeros((n??, 3))
    times = similar(??s)
    directory = "/Users/mattland/Box/work23/20230303-01-hifi_force_for_HSX_coil_1/"
    filename = "20230303-01-hifi_force_for_HSX_coil_1_rtol_$(reltol)_atol_$(abstol)_$(Dates.now()).dat"
    open(directory * filename, "w") do file
        write(file, "??, force_x/length, force_y/length, force_z/length, time\n")
        for j?? in 1:n??
            ?? = ??s[j??]
            println("?? = ", ??)

            #time_data = @timed force = force_finite_thickness(coil, ??; reltol=reltol, abstol=abstol)
            time_data = @timed integrand, force_from_best_fit_circle, force = force_finite_thickness_singularity_subtraction(coil, ??; reltol=reltol, abstol=abstol)
            @show integrand
            @show force_from_best_fit_circle
            @show force
            forces[j??, :] = force
            times[j??] = time_data.time
            println("  time: $(time_data.time)  force: $(force)")
            write(file, "$(??), $(force[1]), $(force[2]), $(force[3]), $(times[j??])\n")
            flush(file)
        end
    end
    println("Finished.")
    
end

function save_high_fidelity_max_B_for_circular_coil_a_scan()
    reltol = 1e-8
    abstol = 1e-8

    a_over_R = 10 .^ collect(((-3.0):(0.0625):(0)))
    println("Values of a/R that will be evaluated: ", a_over_R)

    # Major radius of coil [meters]
    R0 = 1.0
    curve = CurveCircle(R0)

    # Total current [Amperes]
    I = 1.0

    high_fidelity_max_B = similar(a_over_R)
    times = similar(a_over_R)
    for ja in 1:length(a_over_R)
        a = a_over_R[ja] * R0
        println("a = ", a)
        coil = Coil(curve, I, a)

        time_data = @timed modB = hifi_circular_coil_compute_Bz(R0, a, I, R0 - a, 0; reltol=reltol, abstol=abstol)
        high_fidelity_max_B[ja] = modB
        times[ja] = time_data.time
        println("  time: $(time_data.time)  max |B|: $(high_fidelity_max_B[ja])")
    end
    @show high_fidelity_max_B

    directory = "/Users/mattland/Box/work23/20230311-01-circular_coil_max_B_convergence/"
    filename = "circular_coil_max_B_rtol_$(reltol)_atol_$(abstol)_$(Dates.now()).dat"
    open(directory * filename, "w") do file
        write(file, "a/R, max |B|, time\n")
        for ja in 1:length(a_over_R)
            write(file, "$(a_over_R[ja]), $(high_fidelity_max_B[ja]), $(times[ja])\n")
        end
    end
    
end

