"""
For a circular filament coil, show that the IxB force diverges logarithmically
if the evaluation point is on the coil and you merely skip the singular point.
"""
function plot_force_non_convergence_skipping_point_circular()
    # Major radius of coil [meters]
    R0 = 3.0

    # Minor radius of coil [meters]
    aminor = 0.01

    # Total current [Amperes]
    I = 1.0e6

    curve = CurveCircle(R0)
    coil = CoilCircularXSection(curve, I, aminor)
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

function plot_force_non_convergence_skipping_point_circular_with_ours()
    # Major radius of coil [meters]
    R0 = 3.0

    # Minor radius of coil [meters]
    aminor = 0.01

    # Total current [Amperes]
    I = 1.0e6

    curve = CurveCircle(R0)
    coil = CoilCircularXSection(curve, I, aminor)
    r_eval = [R0, 0, 0]

    # Generate numbers of quadrature points to try:
    nns = 25
    ns = [Int(round(10 ^ x)) for x in range(0.5, 4.0, length=nns)]

    force_per_unit_length_skipping_point = zeros(nns)
    force_per_unit_length_ours = zeros(nns)
    for jn in 1:nns
        B = B_filament_fixed(coil, r_eval, ns[jn], drop_first_point=true)
        force_per_unit_length_skipping_point[jn] = I * B[3]
        B = B_singularity_subtraction_fixed(coil, 0, ns[jn])
        force_per_unit_length_ours[jn] = I * B[3]
    end
    @show force_per_unit_length_ours

    x = [minimum(ns), maximum(ns)]

    linewidth = 2
    #scatter(ns, force_per_unit_length, xscale=:log10)
    scatter(
        ns,
        force_per_unit_length_skipping_point,
        xscale=:log10,
        label="Skipping singular point",
        size=(600, 500),
        leg=false,
        c=:blue,
        ms=3,
        msw=0,
        xtickfont=font(10),
        ytickfont=font(10),
        xguidefontsize=12,
        yguidefontsize=12,
        minorgrid=true,
    )
    #plot!(
    #    ns,
    #    force_per_unit_length_skipping_point,
    #    lw=linewidth,
    #    c=:blue,
    #)
    plot!(x, analytic_force_per_unit_length(coil) * [1, 1], label="Analytic", c=:darkorange, lw=5)
    #plot!(x, interpolated_force_per_unit_length(coil) * [1, 1], linestyle=:dash, label="5D integral", c=:green, lw=3)
    #plot!(x, interpolated_force_per_unit_length(coil) * [1, 1], linestyle=:dash, label="5D integral", c=:green, lw=3)
    plot!(x, interpolated_force_per_unit_length(coil) * [1, 1], label="5D integral", c=:green, lw=2)
    #plot!(ns, force_per_unit_length_ours, linestyle=:dash, label="Filament, our method", c=:red, lw=1)
    scatter!(ns, force_per_unit_length_ours, ms=3, msw=0, label="Filament, our method", c=:red, lw=1)
    xlabel!("Number of grid points for filament calculations")
    ylabel!("Force per unit length [N / m]")
    xticks!(10 .^ (1:4))
    xlims!(3, 1e4)
    ylims!((0, Inf))
    #title!("Merely skipping the singular grid point for a filament         \ngives a non-convergent result with O(1) error        ")
    annotate!(8, 2.47e5, text("Analytic", :darkorange))
    annotate!(4.0e3, 2.23e5, text("5D integral", :green))
    annotate!(30, 2.23e5, text("Our method", :red))
    annotate!(300, 0.9e5, text("Filament, skipping singular point", :blue))
    savefig("circular_coil_force_non_convergence_skipping_point.pdf")
end

"""
For a circular filament coil, show that |B| diverges logarithmically
if the evaluation point is on the coil and you merely skip the singular point.
"""
function plot_modB_non_convergence_skipping_point_circular_with_ours()
    # Major radius of coil [meters]
    R0 = 3.0

    # Minor radius of coil [meters]
    aminor = 0.01

    # Total current [Amperes]
    I = 1.0e6

    curve = CurveCircle(R0)
    coil = CoilCircularXSection(curve, I, aminor)
    r_eval = [R0, 0, 0]

    modB_analytic = (
        μ0 * I / (2π * aminor)
        + μ0 * I / (8π * R0) * (6 * log(2.0) - 0.5 + 2 * log(R0 / aminor))
    )

    @time modB_3d_integral = hifi_circular_coil_compute_Bz(R0, aminor, I, R0 - aminor, 0; reltol=1e-6, abstol=1e-6)

    println("modB_analytic:    ", modB_analytic)
    println("modB_3d_integral: ", modB_3d_integral)
    @show relative_difference = (modB_analytic - modB_3d_integral) / (0.5 * (modB_analytic + modB_3d_integral))
    
    # Generate numbers of quadrature points to try:
    nns = 25
    ns = [Int(round(10 ^ x)) for x in range(0.5, 4.0, length=nns)]

    modB_skipping_point = zeros(nns)
    modB_ours = zeros(nns)
    for jn in 1:nns
        B = B_filament_fixed(coil, r_eval, ns[jn], drop_first_point=true)
        modB_skipping_point[jn] = B[3]
        B_filament = B_singularity_subtraction_fixed(coil, 0, ns[jn])
        # Take the b = e_z component of eq (123)-(125) in 
        # 20230326-01_B_in_conductor_for_a_noncircular_finite_thickness_coil.pdf
        # for θ = 0:
        modB_ours[jn] = (
            μ0 * I / (2π * aminor)
            + μ0 * I / (2π * R0) * (0.75 - 1 + 0.5)
            + B_filament[3]
        )
    end
    @show modB_skipping_point

    x = [minimum(ns), maximum(ns)]

    linewidth = 2
    scatter(
        ns,
        modB_skipping_point,
        xscale=:log10,
        yscale=:log10,
        label="Skipping singular point",
        size=(600, 500),
        leg=false,
        c=:blue,
        ms=3,
        msw=0,
        xtickfont=font(10),
        ytickfont=font(10),
        titlefontsize=14,
        xguidefontsize=12,
        yguidefontsize=12,
        titlefonthalign=:left,
        minorgrid=true,
    )
    plot!(x, modB_analytic * [1, 1], label="Analytic", c=:darkorange, lw=5)
    plot!(x, modB_3d_integral * [1, 1], label="3D integral", c=:green, lw=2)
    scatter!(ns, modB_ours, ms=3, msw=0, label="Filament, our method", c=:red, lw=1)
    xlabel!("Number of grid points for filament calculations")
    ylabel!("maximum |B| in the conductor [T]")
    xticks!(10 .^ (1:4))
    xlims!(3, 1e4)
    ylims!((0.03, 40))
    title!("Merely skipping the singular grid point for a filament         \ngives a non-convergent result with significant error        ")
    annotate!(8, 26, text("Analytic", :darkorange))
    annotate!(2.0e3, 26, text("3D integral", :green))
    annotate!(60, 15.0, text("Our method", :red))
    annotate!(220, 0.07, text("Filament, skipping singular point", :blue))
    savefig("circular_coil_modB_non_convergence_skipping_point.pdf")
end


function plot_force_for_HSX()
    coil_num = 2
    curve = get_curve("hsx", coil_num)

    current = -1.5e5

    # minor radius of conductor:
    a = 0.02

    coil = CoilCircularXSection(curve, current, a)
    regularization = a * a / sqrt(ℯ)

    # Number of grid points on which to report the force.
    # (For computing the force, adaptive quadrature will be used.)
    nϕ = 200

    ϕ = (collect(1:nϕ) .- 1) * 2π / nϕ
    force_per_unit_length = zeros(nϕ, 3)
    for j in 1:nϕ
        position, tangent = position_and_tangent(curve, ϕ[j])
        B = B_filament_adaptive(coil, position, regularization=regularization)
        force_per_unit_length[j, :] = current * cross(tangent, B)
    end

    plot(ϕ, force_per_unit_length[:, 1], label="x")
    plot!(ϕ, force_per_unit_length[:, 2], label="y")
    plot!(ϕ, force_per_unit_length[:, 3], label="x")
    xlabel!("Coil parameter ϕ")
    title!("Force per unit length on HSX coil $(coil_num) [N/m]")
end

function plot_integrand()
    coil_num = 1
    curve = get_curve("hsx", coil_num)
    # curve = CurveCircle(2.2)

    current = -1.5e5

    # minor radius of conductor:
    a = 0.01

    # point at which to evaluate the force:
    ϕ0 = 0.0
    ϕ0 = 2π/5

    coil = CoilCircularXSection(curve, current, a)
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

    data = γ_and_3_derivatives(curve, ϕ0)
    r_eval = @view data[:, 1]
    r_prime = @view data[:, 2]
    r_prime_prime = @view data[:, 3]
    r_prime_prime_prime = @view data[:, 4]

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
        plot!(ϕp, integrand[:, j], label=xyz[j:j])
        plot!(ϕp, smoothed_integrand[:, j], label="smoothed, " * xyz[j:j])
        #plot!(ϕp, smoother_integrand[:, j], label="smoother, " * xyz[j:j], ls=:dot)
        #plot!(ϕp, smoothest_integrand[:, j], label="smoothest, " * xyz[j:j], ls=:dash)
    end
    xlabel!("Coil parameter ϕ")
    title!("Regularized Biot-Savart integrand for HSX coil $(coil_num) at ϕ=$(ϕ0)")

end

function plot_force_convergence_single()
    coil_num = 1

    # point at which to evaluate the force:
    #ϕ0 = 0.0
    ϕ0 = 2π/5

    curve = get_curve("hsx", coil_num)
    # curve = CurveCircle(2.2)

    current = -1.5e5

    # minor radius of conductor:
    a = 0.01

    coil = CoilCircularXSection(curve, current, a)
    δ = a * a / sqrt(ℯ)

    # Generate numbers of quadrature points to try:
    nns = 40
    ns = [Int(round(10 ^ x)) for x in range(1.0, 4.0, length=nns)]

    r_eval, tangent0 = position_and_tangent(curve, ϕ0)
    force_per_unit_length = zeros(nns)
    force_per_unit_length_singularity_subtraction = zeros(nns)
    for jn in 1:nns
        B = B_filament_fixed(coil, r_eval, ns[jn], regularization=δ; ϕ_shift=ϕ0)
        #B = B_filament_fixed(coil, r_eval, ns[jn], regularization=δ)
        force_per_unit_length[jn] = current * norm(cross(tangent0, B))

        B = B_singularity_subtraction_fixed(coil, ϕ0, ns[jn])
        force_per_unit_length_singularity_subtraction[jn] = current * norm(cross(tangent0, B))
    end

    scatter(ns, abs.(force_per_unit_length), xscale=:log10, label="original")
    scatter!(ns, abs.(force_per_unit_length_singularity_subtraction), label="singularity subtraction")
    xlabel!("number of quadrature points")
    title!("Force per unit length [N / m]")
end

function plot_force_convergence_single_for_talk()
    coil_num = 1

    # point at which to evaluate the force:
    ϕ0 = 0.0
    #ϕ0 = 4 * 2π/5

    curve = get_curve("hsx", coil_num)
    # curve = CurveCircle(2.2)

    current = -1.5e5

    # minor radius of conductor:
    a = 0.01

    coil = CoilCircularXSection(curve, current, a)
    δ = a * a / sqrt(ℯ)

    # Generate numbers of quadrature points to try:
    nns = 60
    ns = [Int(round(10 ^ x)) for x in range(1.0, 3.0, length=nns)]

    r_eval, tangent0 = position_and_tangent(curve, ϕ0)
    force_per_unit_length = zeros(nns)
    force_per_unit_length_singularity_subtraction = zeros(nns)
    for jn in 1:nns
        B = B_filament_fixed(coil, r_eval, ns[jn], regularization=δ)
        force_per_unit_length[jn] = current * norm(cross(tangent0, B))
        #force_per_unit_length[jn] = current * cross(tangent0, B)[1]

        B = B_singularity_subtraction_fixed(coil, ϕ0, ns[jn])
        force_per_unit_length_singularity_subtraction[jn] = current * norm(cross(tangent0, B))
        #force_per_unit_length_singularity_subtraction[jn] = current * cross(tangent0, B)[1]
    end

    scalefontsizes()
    scalefontsizes(1.3)

    factor = 1e-3
    marker_size = 3.5
    scatter(
        ns, 
        abs.(force_per_unit_length) * factor, 
        xscale=:log10, 
        ms=marker_size,
        msw=0,
        #label="original",
        label=false,
        framestyle=:box,
        titlefontsize=15,
        minorgrid=true,
    )
    scatter!(
        ns, 
        abs.(force_per_unit_length_singularity_subtraction) * factor, 
        #label="singularity subtraction",
        label=false,
        ms=marker_size,
        msw=0,
    )
    xlabel!("Number of quadrature points")
    title!("Force per unit length |dF/dℓ| [kN / m]")
    ylims!(0, 60)
    xlims!(10, 1e3)
    savefig("/Users/mattland/Box/work23/20230511-01-HSX_force_convergence_vs_nquadpoints.pdf")
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
        coil = CoilCircularXSection(curve, current, a)
    
        for jϕ in eachindex(ϕs)
            ϕ0 = ϕs[jϕ]

            r_eval, tangent0 = position_and_tangent(curve, ϕ0)
            force_per_unit_length = zeros(nns)
            force_per_unit_length_singularity_subtraction = zeros(nns)
            for jn in 1:nns
                B = B_filament_fixed(coil, r_eval, ns[jn], regularization=δ; ϕ_shift=ϕ0)
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
    coil = CoilCircularXSection(curve, I, aminor)

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
    for ja in eachindex(a_over_R)
        a = a_over_R[ja] * R0
        println("a = ", a)
        coil = CoilCircularXSection(curve, I, a)

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
        for ja in eachindex(a_over_R)
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
    for ja in eachindex(a_over_R)
        a = a_over_R[ja]
        coil = CoilCircularXSection(curve, I, a)

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
        for ja in eachindex(a_over_R)
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

    coil = CoilCircularXSection(curve, I, a_over_R)

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

function save_high_fidelity_force_for_circular_coil2_tol_scan()
    a_over_R = 0.01

    # Major radius of coil [meters]
    R0 = 1.0
    curve = CurveCircle(R0)

    # Total current [Amperes]
    I = 1.0
    
    coil = CoilCircularXSection(curve, I, a_over_R)

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
    coil = CoilCircularXSection(curve, current, aminor)
    regularization = aminor * aminor / sqrt(exp(1))

    force = zeros((neval, 3))
    force_circular_approximation = zeros((neval, 3))
    ϕ = [2π * (j - 1) / neval for j in 1:neval]
    for j in 1:neval
        position, tangent_vector = position_and_tangent(curve, ϕ[j])
        B = B_filament_adaptive(coil, position; regularization=regularization)
        force[j, :] = current * cross(tangent_vector, B)
        force_circular_approximation[j, :] = force_locally_circular_approximation(coil, ϕ[j])
    end

    plot(ϕ, force[:, 1], label="filament, adaptive quadrature")
    plot!(ϕ, force_circular_approximation[:, 1], label="locally circular approximation")
    #plot!(ϕ, force[:, 2], label="y")
    #plot!(ϕ, force[:, 3], label="z")
    xlabel!("ϕ")
    title!("dFₓ/dℓ for HSX coil 1")
end

function plot_Fx_for_HSX_coil_1()
    # Number of evaluation points for the force:
    neval = 300

    curve = get_curve("hsx", 1)
    current = 1.0
    # Data from HsxCoilsNescToFinite.pdf:
    x_sectional_area = (56.8e-3) * (129.6e-3)
    aminor = sqrt(x_sectional_area / π)
    @show aminor

    coil = CoilCircularXSection(curve, current, aminor)
    regularization = aminor * aminor / sqrt(exp(1))

    n_quadpoints_arr = [100, 1000, 10000]
    n_quadpointss = length(n_quadpoints_arr)
    force = zeros((neval, 3))
    force_circular_approximation = zeros((neval, 3))
    force_drop_singular_point = zeros((neval, 3, n_quadpointss))
    ϕ = [2π * (j - 1) / neval for j in 1:neval]
    for j in 1:neval
        position, tangent_vector = position_and_tangent(curve, ϕ[j])
        B = B_filament_adaptive(coil, position; regularization=regularization)
        force[j, :] = current * cross(tangent_vector, B)
        force_circular_approximation[j, :] = force_locally_circular_approximation(coil, ϕ[j])
        for jquadpoints in 1:n_quadpointss
            n_quadpoints = n_quadpoints_arr[jquadpoints]
            dϕ = 2π / n_quadpoints
            B = [0.0, 0.0, 0.0]
            for k in 2:n_quadpoints
                ϕ_shifted = ϕ[j] + (k - 1) * dϕ
                #B += d_B_d_ϕ(coil, ϕ_shifted, position, regularization=regularization)
                B += d_B_d_ϕ(coil, ϕ_shifted, position, regularization=0)
            end
            B *= dϕ
            force_drop_singular_point[j, :, jquadpoints] = current * cross(tangent_vector, B)
        end
    end

    plot(ϕ, force[:, 1], label="filament, adaptive quadrature")
    #plot!(ϕ, force_circular_approximation[:, 1], label="locally circular approximation")
    for jquadpoints in 1:n_quadpointss
        plot!(ϕ, force_drop_singular_point[:, 1, jquadpoints], label="filament, drop singular point, $(n_quadpoints_arr[jquadpoints]) points")
    end
    xlabel!("ϕ")
    title!("dFₓ/dℓ for HSX coil 1")
end


function save_high_fidelity_force_for_HSX(;
    aminor = nothing,
    nϕ = 4,
    reltol = 1e-1,
    abstol = 1e-1,
)
    ϕs = [2π * (j - 1) / nϕ for j in 1:nϕ]
    println("Values of ϕ that will be evaluated: ", ϕs)

    curve = get_curve("hsx", 1)
    #curve = CurveCircle(0.5)

    # Total current [Amperes]
    current = 150.0e3
    # Data from HsxCoilsNescToFinite.pdf:
    x_sectional_area = (56.8e-3) * (129.6e-3)
    if aminor === nothing
        aminor = sqrt(x_sectional_area / π)
    end
    coil = CoilCircularXSection(curve, current, aminor)
    @show aminor

    forces = zeros((nϕ, 3))
    forces_with_Siena_trick = zeros((nϕ, 3))
    times = similar(ϕs)
    times_with_Siena_trick = similar(ϕs)
    directory = "/Users/mattland/Box/work23/20230429-01_hifi_force_for_HSX_coil_1/"
    datestr = replace("$(Dates.now())", ":" => ".")
    filename = "20230429-01-hifi_force_for_HSX_coil_1_a$(aminor)_nphi$(nϕ)_rtol$(reltol)_atol$(abstol)_$(datestr).dat"
    open(directory * filename, "w") do file
        write(file, "ϕ, force_x/length, force_y/length, force_z/length, time, force_x/length, force_y/length, force_z/length, time\n")
        for jϕ in 1:nϕ
            ϕ = ϕs[jϕ]
            println("ϕ = ", ϕ)

            println("  Without Siena's trick:")
            time_data = @timed force = force_finite_thickness(coil, ϕ; reltol=reltol, abstol=abstol)
            @show force
            @show time_data.time
            println("  With Siena's trick:")
            time_data_with_Siena_trick = @timed integrand, force_from_best_fit_circle, force_with_Siena_trick = force_finite_thickness_singularity_subtraction(coil, ϕ; reltol=reltol, abstol=abstol)
            @show integrand
            @show force_from_best_fit_circle
            @show force_with_Siena_trick
            @show time_data_with_Siena_trick.time
            forces[jϕ, :] = force
            times[jϕ] = time_data.time
            forces_with_Siena_trick[jϕ, :] = force_with_Siena_trick
            times_with_Siena_trick[jϕ] = time_data_with_Siena_trick.time
            #println("  time: $(time_data.time)  force: $(force)")
            write(file, "$(ϕ), $(force[1]), $(force[2]), $(force[3]), $(times[jϕ]), $(force_with_Siena_trick[1]), $(force_with_Siena_trick[2]), $(force_with_Siena_trick[3]), $(times_with_Siena_trick[jϕ])\n")
            flush(file)
        end
    end
    println("Finished.")
    
end

function save_high_fidelity_force_for_HSX_parallel(;
    aminor = nothing,
    current = 150.0e3,
    nϕ = 4,
    reltol = 1e-1,
    abstol = 1e-1,
)
    ϕs = [2π * (j - 1) / nϕ for j in 1:nϕ]
    println("Values of ϕ that will be evaluated: ", ϕs)

    # Data from HsxCoilsNescToFinite.pdf:
    x_sectional_area = (56.8e-3) * (129.6e-3)
    if aminor === nothing
        aminor = sqrt(x_sectional_area / π)
    end
    @show aminor

    forces = zeros((nϕ, 3))
    times = similar(ϕs)
    println("Number of threads: ", Threads.nthreads())
    directory = "/Users/mattland/Box/work23/20230429-01_hifi_force_for_HSX_coil_1/"
    datestr = replace("$(Dates.now())", ":" => ".")
    filename = "20230429-01-hifi_force_for_HSX_coil_1_a$(aminor)_nphi$(nϕ)_rtol$(reltol)_atol$(abstol)_$(datestr).dat"
    Threads.@threads for jϕ in 1:nϕ
        # CurveXYZFourier has buffers that need to have distinct contents in
        # different threads, so we define the curve here inside the loop, where
        # all new variables are distinct for each thread.
        curve = get_curve("hsx", 1)
        coil = CoilCircularXSection(curve, current, aminor)
        ϕ = ϕs[jϕ]
        println("Thread $(Threads.threadid()) is handling jϕ = $(jϕ) of $(nϕ): ϕ = ", ϕ)
        
        time_data = @timed force = force_finite_thickness(coil, ϕ; reltol=reltol, abstol=abstol)
        forces[jϕ, :] = force
        times[jϕ] = time_data.time
    end
    @show forces
    @show times
    open(directory * filename, "w") do file
        write(file, "ϕ, force_x/length, force_y/length, force_z/length, time\n")
        for jϕ in 1:nϕ
            write(file, "$(ϕs[jϕ]), $(forces[jϕ, 1]), $(forces[jϕ, 2]), $(forces[jϕ, 3]), $(times[jϕ])\n")
        end
    end
    println("Finished.")
    
end

function plot_high_fidelity_force_for_HSX_coil(filename)
    directory = "/Users/mattland/Box/work23/20230429-01_hifi_force_for_HSX_coil_1/"
    f = CSV.File(directory * filename)
    ϕ = f.ϕ

    p1 = plot(ϕ, f[" force_x/length"], label="direct")
    plot!(ϕ, f[" force_x/length_1"], linestyle=:dash, label="Siena")
    xlabel!("ϕ")
    ylabel!("force_x/length")

    p2 = plot(ϕ, f[" force_y/length"], label="direct")
    plot!(ϕ, f[" force_y/length_1"], linestyle=:dash, label="Siena")
    xlabel!("ϕ")
    ylabel!("force_y/length")
    
    p3 = plot(ϕ, f[" force_z/length"], label="direct")
    plot!(ϕ, f[" force_z/length_1"], linestyle=:dash, label="Siena")
    xlabel!("ϕ")
    ylabel!("force_z/length")
    
    p4 = plot(ϕ, f[" time"], label="direct")
    plot!(ϕ, f[" time_1"], linestyle=:dash, label="Siena")
    xlabel!("ϕ")
    ylabel!("time")
    
    plot(p1, p2, p3, p4, layout=(2, 2), dpi=100, size=(700, 600))

end

function plot_high_fidelity_force_for_HSX_coil_compare(filenames)
    directory = "/Users/mattland/Box/work23/20230429-01_hifi_force_for_HSX_coil_1/"
    fs = [CSV.File(directory * filename) for filename in filenames]
    n = length(fs)

    scalefontsizes()
    scalefontsizes(0.5)

    p1 = plot()
    for j in 1:n 
        plot!(fs[j].ϕ, fs[j][" force_x/length"], label=false)
    end
    xlabel!("ϕ")
    ylabel!("force_x/length")

    p2 = plot()
    for j in 1:n 
        plot!(fs[j].ϕ, fs[j][" force_y/length"], label=false)
    end
    xlabel!("ϕ")
    ylabel!("force_y/length")
    
    p3 = plot()
    for j in 1:n 
        plot!(fs[j].ϕ, fs[j][" force_z/length"], label=false)
    end
    xlabel!("ϕ")
    ylabel!("force_z/length")
    
    p4 = plot()
    for j in 1:n 
        plot!(fs[j].ϕ, fs[j][" time"], label=filenames[j], yscale=:log10)
    end
    xlabel!("ϕ")
    ylabel!("time")
    
    plot(p1, p2, p3, p4, layout=(2, 2), dpi=100, size=(700, 600))

end

function save_high_fidelity_force_for_circular_coil_rectangular_xsection()
    #reltol = 3e-7
    #abstol = 3e-7
    reltol = 1e-3
    abstol = 1e-3

    # geometric mean of a and b [meters]
    d = 0.01

    #b_over_as = [1.0]
    b_over_as = 10 .^ collect(((-1.5):(0.125):(1.5)))
    println("Values of b/a that will be evaluated: ", b_over_as)

    # Major radius of coil [meters]
    R0 = 1.0

    # Total current [Amperes]
    I = 1.0

    as = similar(b_over_as)
    bs = similar(b_over_as)
    force_hifi = similar(b_over_as)
    force_analytic = similar(b_over_as)
    times = similar(b_over_as)
    println("Number of threads: $(Threads.nthreads())")
    Threads.@threads for ja in eachindex(b_over_as)
        b_over_a = b_over_as[ja]
        a = d / sqrt(b_over_a)
        b = d * sqrt(b_over_a)
        as[ja] = a
        bs[ja] = b
        println("a = ", a)
        curve = CurveCircle(R0)
        coil = CoilRectangularXSection(curve, I, a, b, FrameCircle())

        time_data = @timed force = force_finite_thickness(coil, 0.0; reltol=reltol, abstol=abstol)
        force_hifi[ja] = force[1]
        force_analytic[ja] = analytic_force_per_unit_length(coil)
        times[ja] = time_data.time
        println("  thread: $(Threads.threadid())  a: $(a)  b: $(b)  time: $(time_data.time)  analytic force: $(force_analytic[ja])  hifi force: $(force_hifi[ja])  ratio: $(force_analytic[ja] / force_hifi[ja])")
    end

    directory = "/Users/mattland/Box/work23/20230625-01-rectangular_xsection_force/"
    datestr = replace("$(Dates.now())", ":" => ".")
    filename = "circular_coil_rectangular_xsection_force_d$(d)_rtol$(reltol)_atol$(abstol)_$(datestr).dat"
    open(directory * filename, "w") do file
        write(file, "R0, current, d, reltol, abstol\n")
        write(file, "$(R0), $(I), $(d), $(reltol), $(abstol)\n")
        write(file, "a,b,Fx_hifi,Fx_analytic,time\n")
        for ja in eachindex(b_over_as)
            write(file, "$(as[ja]),$(bs[ja]),$(force_hifi[ja]),$(force_analytic[ja]),$(times[ja])\n")
        end
    end
    
end

function plot_high_fidelity_force_for_circular_coil_rectangular_xsection()
    directory = "/Users/mattland/Box/work23/20230625-01-rectangular_xsection_force/"
    filenames = [
        #"circular_coil_rectangular_xsection_force_d0.001_rtol0.01_atol0.01_2023-06-29T13.47.47.338.dat",
        "circular_coil_rectangular_xsection_force_d0.001_rtol0.001_atol0.001_2023-06-29T13.50.25.670.dat",
    ]
    filenames = [
        "circular_coil_rectangular_xsection_force_d0.01_rtol0.01_atol0.01_2023-06-27T20.47.20.198.dat",
        "circular_coil_rectangular_xsection_force_d0.01_rtol0.001_atol0.001_2023-06-27T20.50.34.491.dat",
    ]
    plot()
    common_R0 = -1.0
    common_I = -1.0
    common_d = -1.0
    for j in eachindex(filenames)
        filename = filenames[j]
        file = open(directory * filename, "r")

        # Read and parse header:
        line = readline(file)
        line = readline(file)
        splitline = split(line, ",")
        R0 = parse(Float64, splitline[1])
        I = parse(Float64, splitline[2])
        d = parse(Float64, splitline[3])
        reltol = parse(Float64, splitline[4])
        abstol = parse(Float64, splitline[5])
        close(file)
        # Make sure all of the files being plotted differ only in tolerances:
        if common_R0 < 0
            common_R0 = R0
        else
            @assert common_R0 ≈ R0
        end 
        if common_I < 0
            common_I = I
        else
            @assert common_I ≈ I
        end 
        if common_d < 0
            common_d = d
        else
            @assert common_d ≈ d
        end 

        f = CSV.File(directory * filename, header=3)
        plot!(f.b ./ f.a, f.Fx_hifi, 
            label="hifi tol $(reltol)",
            xscale=:log10
        )
        plot!(f.b ./ f.a, f.Fx_analytic, label="1D tol $(reltol)")
    end
    xlabel!("b / a")
    ylabel!("d F_x / d ℓ")
end

function save_high_fidelity_Bz_for_circular_coil_a_scan()
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
    for ja in eachindex(a_over_R)
        a = a_over_R[ja] * R0
        println("a = ", a)
        coil = CoilCircularXSection(curve, I, a)

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
        for ja in eachindex(a_over_R)
            write(file, "$(a_over_R[ja]), $(high_fidelity_max_B[ja]), $(times[ja])\n")
        end
    end
    
end

function save_high_fidelity_Bz_for_circular_coil_2D()
    reltol = 1e-5
    abstol = 1e-5

    # Resolution for evaluating B
    nz = 58
    nx = 59

    # Major radius of coil [meters]
    R0 = 1.0

    # Min radius of coil [meters]
    a = 1e-2

    # Total current [Amperes]
    I = 1.0

    high_fidelity_Bz = zeros(nz, nx)
    for jz in 1:nz
        z = ((jz - 1.0) / (nz - 1) - 0.5) * 2 * a
        println("z: $(z)")
        for jx = 1:nx
            x = R0 + ((jx - 1.0) / (nx - 1) - 0.5) * 2 * a
            Bz = hifi_circular_coil_compute_Bz(R0, a, I, x, z; reltol=reltol, abstol=abstol)
            high_fidelity_Bz[jz, jx] = Bz
            println("  x: $(x)  Bz: $(Bz)")
        end
    end
    @show high_fidelity_Bz

    directory = "/Users/mattland/Box/work23/20230317-01-circular_coil_Bz_2d/"
    datestr = replace("$(Dates.now())", ":" => ".")
    filename = "circular_coil_Bz_a_$(a)_rtol_$(reltol)_atol_$(abstol)_$(datestr).dat"
    open(directory * filename, "w") do file
        for jz in 1:nz
            write(file, "$(high_fidelity_Bz[jz, 1])")
            for jx = 2:nx
                write(file, ", $(high_fidelity_Bz[jz, jx])")
            end
            write(file, "\n")
        end
    end
    
end

function save_high_fidelity_B_vector_for_circular_coil_2D()
    reltol = 1e-5
    abstol = 1e-5

    # Resolution for evaluating B
    nz = 58
    nx = 59

    # Major radius of coil [meters]
    R0 = 1.0

    # Min radius of coil [meters]
    a = 1e-3

    # Total current [Amperes]
    I = 1.0

    curve = CurveCircle(R0)
    coil = CoilCircularXSection(curve, I, a)
    high_fidelity_B = zeros(nz, nx, 3)
    for jz in 1:nz
        z = ((jz - 1.0) / (nz - 1) - 0.5) * 2 * a
        println("z: $(z)")
        for jx = 1:nx
            x = R0 + ((jx - 1.0) / (nx - 1) - 0.5) * 2 * a
            B = B_finite_thickness(coil, [x, 0, z]; reltol=reltol, abstol=abstol)
            high_fidelity_B[jz, jx, :] = B
            println("  x: $(x)  B: $(B)")
        end
    end
    @show high_fidelity_B

    directory = "/Users/mattland/Box/work23/20230317-01-circular_coil_Bz_2d/"
    filename = "circular_coil_B_vector_a_$(a)_rtol_$(reltol)_atol_$(abstol)_$(Dates.now()).dat"
    open(directory * filename, "w") do file
        for jxyz in 1:3
            for jz in 1:nz
                write(file, "$(high_fidelity_B[jz, 1, jxyz])")
                for jx = 2:nx
                    write(file, ", $(high_fidelity_B[jz, jx, jxyz])")
                end
                write(file, "\n")
            end
        end
    end
    
end

function save_high_fidelity_B_vector_for_circular_coil_rectangular_xsection()
    reltol = 1e-5
    abstol = 1e-5

    # Resolution for evaluating B
    nz = 25
    nx = 26

    # Major radius of coil [meters]
    R0 = 1.0

    # Cross-section dimensions [meters]
    a = 2e-3
    b = 1e-3

    # Total current [Amperes]
    I = 1.0

    curve = CurveCircle(R0)
    coil = CoilRectangularXSection(curve, I, a, b, FrameCircle())
    high_fidelity_B = zeros(nz, nx, 3)
    for jz in 1:nz
        z = ((jz - 1.0) / (nz - 1) - 0.5) * b
        println("z: $(z)")
        for jx = 1:nx
            x = R0 + ((jx - 1.0) / (nx - 1) - 0.5) * a
            B = B_finite_thickness(coil, [x, 0, z]; reltol=reltol, abstol=abstol)
            high_fidelity_B[jz, jx, :] = B
            println("  x: $(x)  B: $(B)")
        end
    end
    @show high_fidelity_B

    directory = "/Users/mattland/Box/work23/20230621-01-rectangular_xsection_B/"
    datestr = replace("$(Dates.now())", ":" => ".")
    filename = "circular_coil_rectangular_xsection_B_vector_a$(a)_b$(b)_nx$(nx)_nz$(nz)_rtol$(reltol)_atol$(abstol)_$(datestr).dat"
    open(directory * filename, "w") do file
        write(file, "R0, a, b, I, reltol, abstol, nx, nz\n")
        write(file, "$(R0), $(a), $(b), $(I), $(reltol), $(abstol), $(nx), $(nz)\n")
        for jxyz in 1:3
            for jz in 1:nz
                for jx = 1:nx
                    write(file, "$(high_fidelity_B[jz, jx, jxyz])\n")
                end
            end
        end
    end
    
end

function plot_B_vector_for_circular_coil_rectangular_xsection(
    filename = "circular_coil_rectangular_xsection_B_vector_a0.001_b0.002_nx26_nz25_rtol1.0e-5_atol1.0e-5_2023-06-21T13.35.24.678.dat"
)
    directory = "/Users/mattland/Box/work23/20230621-01-rectangular_xsection_B/"
    file = open(directory * filename, "r")

    # Read and parse header:
    line = readline(file)
    line = readline(file)
    splitline = split(line, ",")
    R0 = parse(Float64, splitline[1])
    a = parse(Float64, splitline[2])
    b = parse(Float64, splitline[3])
    I = parse(Float64, splitline[4])
    reltol = parse(Float64, splitline[5])
    abstol = parse(Float64, splitline[6])
    nx = parse(Int, splitline[7])
    nz = parse(Int, splitline[8])
    println("Read nx=$(nx), nz=$(nz)")

    curve = CurveCircle(R0)
    coil = CoilRectangularXSection(curve, I, a, b, FrameCircle())

    # Now read the main B data:
    high_fidelity_B = zeros(nz, nx, 3)
    for jxyz in 1:3
        for jz in 1:nz
            for jx = 1:nx
                high_fidelity_B[jz, jx, jxyz] = parse(Float64, readline(file))
            end
        end
    end
    close(file)
    #@show high_fidelity_B

    # Set up grid of subplots
    n_plots = 3 * 3
    n_cols = Int(ceil(0.9 * sqrt(n_plots)))
    n_rows = Int(ceil(n_plots / n_cols))
    @show n_plots, n_rows, n_cols
    xyz = "xyz"

    layout = (n_rows, n_cols)
    plots = Array{Any, 1}(undef, n_plots)
    scalefontsizes()
    scalefontsizes(0.5)

    Plots.gr_cbar_width[] = 0.005

    almost_one = 1 - (1e-6)
    u = [-((jx - 1) / (nx - 1) * 2 - 1) * almost_one for jx in 1:nx]
    v = [((jz - 1) / (nz - 1) * 2 - 1) * almost_one for jz in 1:nz]
    u2d = [u[jx] for jz in 1:nz, jx in 1:nx]
    v2d = [v[jz] for jz in 1:nz, jx in 1:nx]
    x = R0 .- u * a / 2
    z = v * b / 2
    κ1 = 1 / R0
    κ2 = 0

    B0 = zeros(nz, nx, 3)
    for jx in 1:nx
        for jz in 1:nz
            B0_p, B0_q = CoilForces.rectangular_xsection_B0(coil, u[jx], v[jz])
            B0[jz, jx, 1] = -B0_p
            B0[jz, jx, 3] = B0_q
        end
    end

    position = [R0, 0.0, 0.0]
    B_regularized_plus_extra_term = B_filament_adaptive(coil, position; regularization=compute_regularization(coil))
    B_regularized_plus_extra_term[3] += μ0 * coil.current / (8π * R0) * (4 + 2 * log(2) + log(CoilForces.rectangular_xsection_δ(a, b)))
    for subtract_leading_order in [false, true]
        index = 1
        
        for jxyz in 1:3
            B_analytic = zeros(nz, nx)
            for jz in 1:nz
                for jx in 1:nx
                    Bκ_p, Bκ_q = CoilForces.rectangular_xsection_Bκ(coil, κ1, κ2, u[jx], v[jz])
                    Bκ_term = [-Bκ_p, 0, Bκ_q]
                    B_analytic[jz, jx] = (
                        B_regularized_plus_extra_term[jxyz] 
                        + B0[jz, jx, jxyz]
                        + Bκ_term[jxyz]
                        )
                end
            end
            B_hifi = high_fidelity_B[:, :, jxyz]
            if subtract_leading_order
                B_analytic -= B0[:, :, jxyz]
                B_hifi -= B0[:, :, jxyz]
            end
            for variant in ["hifi", "analytic", "difference"]
                if variant == "hifi"
                    data = B_hifi
                elseif variant == "analytic"
                    data = B_analytic
                elseif variant == "difference"
                    data = B_hifi - B_analytic
                end

                #maxB = maximum(leading_order_solution)
                #minB = minimum(leading_order_solution)
                #if maxB == minB
                #    maxB += 1e-10
                #end
                #contour_levels = collect(range(minB, maxB, length=25))
                #@show contour_levels

                @show 
                plots[index] = contour(x, z, data,
                    aspect_ratio = :equal,
                    #levels=contour_levels,
                )
                index += 1
                title_str = @sprintf "B%s [Tesla]" xyz[jxyz]
                #title!("B$(xyz[jxyz]) [Tesla] at ϕ=$(ϕ)")

                #if hifi
                #    title_str = "HiFi " * title_str
                #else
                #    title_str = "analytic " * title_str
                #end
                title_str = variant * " " * title_str
                title!(title_str)
                xlabel!("x [meters]")
                ylabel!("z [meters]")
            end
        end
        plot(plots..., layout=layout, dpi=100, size=(1100, 850))
        if subtract_leading_order
            filename_extension = "_without_leading_order"
        else
            filename_extension = ""
        end
        #annotate!((0., 0., directory * filename), subplot=1)
        savefig(directory * "circular_coil_rectangular_xsection_B_vector" * filename_extension * ".pdf")
    end

end

function save_high_fidelity_B_vector_for_HSX_coil(;
    aminor = 0.01,  # Minor radius of coil [meters]
    I = 1.0e4,  # Total current [Amperes]
    reltol = 1e-4,
    abstol = 1e-4,
    nϕ = 4,
    nn = 2,
    nb = 3,
)
    curve = get_curve("hsx", 1)
    #curve = CurveCircle(1.0)
    coil = CoilCircularXSection(curve, I, aminor)
    high_fidelity_B = zeros(nϕ, nn, nb, 3)
    for jϕ in 1:nϕ
        println("Processing jϕ = $(jϕ) of $(nϕ)")
        ϕ = 2π * (jϕ - 1) / nϕ
        differential_arclength, curvature, torsion, position, tangent, normal, binormal = Frenet_frame(curve, ϕ)
        for jn in 1:nn
            println("  Processing jn = $(jn) of $(nn)")
            xn = ((jn - 1.0) / (nn - 1) - 0.5) * 2 * aminor
            for jb = 1:nb
                xb = ((jb - 1.0) / (nb - 1) - 0.5) * 2 * aminor
                eval_point = position + xn * normal + xb * binormal
                θ_shift = atan(xb, xn)
                #θ_shift = atan(xn, xb)
                #println("      About to evaluate at ", eval_point, " ϕ=", ϕ, " θ_shift=", θ_shift)
                B = B_finite_thickness(coil, eval_point; reltol=reltol, abstol=abstol, ϕ_shift=ϕ, θ_shift=θ_shift)
                #B = B_finite_thickness(coil, eval_point; reltol=reltol, abstol=abstol)
                high_fidelity_B[jϕ, jn, jb, :] = B
                println("    jb: $(jb)  B: $(B)")
            end
        end
    end
    @show high_fidelity_B
    #return

    directory = "/Users/mattland/Box/work23/20230426-01-HSX_coil_hifi_B_vector/"
    datestr = replace("$(Dates.now())", ":" => ".")
    filename = "HSX_B_a$(aminor)_rtol$(reltol)_atol$(abstol)_nphi$(nϕ)_nn$(nn)_nb$(nb)_$(datestr).dat"
    open(directory * filename, "w") do file
        write(file, "aminor, I, reltol, abstol, nϕ, nn, nb\n")
        write(file, "$(aminor), $(I), $(reltol), $(abstol), $(nϕ), $(nn), $(nb)\n")
        for jϕ in 1:nϕ
            for jn in 1:nn
                for jb = 1:nb
                    for jxyz in 1:3
                        write(file, "$(high_fidelity_B[jϕ, jn, jb, jxyz])\n")
                    end
                end
            end
        end
    end
    
end

function save_high_fidelity_B_vector_for_HSX_coil_rectangular_xsection(;
    a = 0.1296,
    b = 0.0568,
    I = 150.0e3,  # Total current [Amperes]
    reltol = 1e-3,
    abstol = 1e-3,
    nϕ = 4,
    nu = 2,
    nv = 3,
)
    high_fidelity_B = zeros(nϕ, nu, nv, 3)
    println("Number of threads: ", Threads.nthreads())
    Threads.@threads for ju in 1:nu
        curve = get_curve("hsx", 1)
        #curve = CurveCircle(1.0)
        coil = CoilRectangularXSection(curve, I, a, b, FrameCentroid(curve))
        println("Thread $(Threads.threadid()) is processing ju = $(ju) of $(nu)")
        u = ((ju - 1.0) / (nu - 1) - 0.5) * 2
        for jϕ in 1:nϕ
            println("  Processing jϕ = $(jϕ) of $(nϕ)")
            ϕ = 2π * (jϕ - 1) / nϕ
            differential_arclength, curvature, torsion, position, tangent, normal, binormal = Frenet_frame(curve, ϕ)
            p, q = get_frame(coil.frame, ϕ, position, tangent, normal)
                for jv = 1:nv
                v = ((jv - 1.0) / (nv - 1) - 0.5) * 2
                eval_point = position + (u * a / 2) * p + (v * b / 2) * q
                #println("      About to evaluate at ", eval_point, " ϕ=", ϕ, " θ_shift=", θ_shift)
                B = B_finite_thickness(coil, eval_point; reltol=reltol, abstol=abstol, ϕ_shift=ϕ)
                high_fidelity_B[jϕ, ju, jv, :] = B
                println("    jv: $(jv)  B: $(B)")
            end
        end
    end
    @show high_fidelity_B
    #return

    directory = "/Users/mattland/Box/work23/20230621-01-rectangular_xsection_B/"
    datestr = replace("$(Dates.now())", ":" => ".")
    filename = "HSX_rectangular_xsection_B_a$(a)_b$(b)_rtol$(reltol)_atol$(abstol)_nphi$(nϕ)_nu$(nu)_nv$(nv)_$(datestr).dat"
    open(directory * filename, "w") do file
        write(file, "a, b, I, reltol, abstol, nϕ, nu, nv\n")
        write(file, "$(a), $(b), $(I), $(reltol), $(abstol), $(nϕ), $(nu), $(nv)\n")
        for jϕ in 1:nϕ
            for ju in 1:nu
                for jv = 1:nv
                    for jxyz in 1:3
                        write(file, "$(high_fidelity_B[jϕ, ju, jv, jxyz])\n")
                    end
                end
            end
        end
    end
    
end

function plot_high_fidelity_B_vector_for_HSX_coil(
    filename="HSX_B_a0.01_rtol1.0e-5_atol1.0e-5_nphi4_nn2_nb3_2023-04-26T13.25.33.056.dat"
)
    directory = "/Users/mattland/Box/work23/20230426-01-HSX_coil_hifi_B_vector/"
    file = open(directory * filename, "r")

    # Read and parse header:
    line = readline(file)
    line = readline(file)
    splitline = split(line, ",")
    aminor = parse(Float64, splitline[1])
    I = parse(Float64, splitline[2])
    nϕ = parse(Int, splitline[5])
    nn = parse(Int, splitline[6])
    nb = parse(Int, splitline[7])
    println("Read nϕ=$(nϕ), nn=$(nn), nb=$(nb)")

    curve = get_curve("hsx", 1)
    #curve = CurveCircle(1.0)
    coil = CoilCircularXSection(curve, I, aminor)
    regularization = aminor * aminor / sqrt(ℯ)

    # Now read the main B data:
    high_fidelity_B = zeros(nϕ, nn, nb, 3)
    for jϕ in 1:nϕ
        for jn in 1:nn
            for jb = 1:nb
                for jxyz in 1:3
                    high_fidelity_B[jϕ, jn, jb, jxyz] = parse(Float64, readline(file))
                end
            end
        end
    end
    close(file)
    #@show high_fidelity_B

    # Set up grid of subplots
    n_plots = nϕ * 3 * 2
    n_cols = Int(ceil(0.9 * sqrt(n_plots)))
    n_rows = Int(ceil(n_plots / n_cols))
    @show n_plots, n_rows, n_cols
    xyz = "xyz"

    layout = (n_rows, n_cols)
    plots = Array{Any, 1}(undef, n_plots)
    scalefontsizes()
    scalefontsizes(0.5)

    Plots.gr_cbar_width[] = 0.005

    uplot = [((jn - 1) / (nn - 1) * 2 - 1) * aminor for jn in 1:nn]
    vplot = [((jb - 1) / (nb - 1) * 2 - 1) * aminor for jb in 1:nb]
    u2d = [uplot[jn] for jb in 1:nb, jn in 1:nn]
    v2d = [vplot[jb] for jb in 1:nb, jn in 1:nn]
    ρ = @. sqrt(u2d^2 + v2d^2) / aminor
    θ = atan.(v2d, u2d)
    @show ρ
    @show θ

    for subtract_leading_order in [false, true]
        index = 1
        for jϕ in 1:nϕ
            ϕ = 2π * (jϕ - 1) / nϕ
            differential_arclength, curvature, torsion, position, tangent, normal, binormal = Frenet_frame(curve, ϕ)

            B_regularized = B_filament_adaptive(coil, position; regularization=regularization)
            @show B_regularized
            
            for jxyz in 1:3
                leading_order_solution = μ0 * I / (2π * aminor * aminor) * (
                    -normal[jxyz] * v2d + binormal[jxyz] * u2d
                )
                for hifi in [true, false]
                    if hifi
                        data = high_fidelity_B[jϕ, :, :, jxyz]'
                    else
                        data = zeros(nb, nn)
                        for jb in 1:nb
                            for jn in 1:nn
                                data[jb, jn] = (
                                    B_regularized[jxyz] 
                                    + CoilForces.B_local(coil, curvature, normal[jxyz], binormal[jxyz], ρ[jb, jn], θ[jb, jn])
                                )
                            end
                        end
                    end
                    if subtract_leading_order
                        data -= leading_order_solution
                    end

                    # Don't plot data outside the coil - make those points NaN
                    for jb in 1:nb
                        for jn in 1:nn
                            if uplot[jn]^2 + vplot[jb]^2 > aminor^2
                                data[jb, jn] = NaN
                            end
                        end
                    end

                    maxB = maximum(leading_order_solution)
                    minB = minimum(leading_order_solution)
                    if maxB == minB
                        maxB += 1e-10
                    end
                    contour_levels = collect(range(minB, maxB, length=25))
                    #@show contour_levels

                    plots[index] = contour(uplot, vplot, data,
                        aspect_ratio = :equal,
                        #levels=contour_levels,
                    )
                    index += 1
                    title_str = @sprintf "B%s [Tesla] at ϕ=%.2f" xyz[jxyz] ϕ
                    #title!("B$(xyz[jxyz]) [Tesla] at ϕ=$(ϕ)")
                    if hifi
                        title_str = "HiFi " * title_str
                    else
                        title_str = "analytic " * title_str
                    end
                    title!(title_str)
                    xlabel!("u [meters]")
                    ylabel!("v [meters]")
                    nθ = 150
                    θplot = collect(range(0, 2π, length=nθ))
                    plot!(aminor * cos.(θplot), aminor * sin.(θplot), linewidth=1.5, color=:black, label=nothing)
                end
            end
        end
        plot(plots..., layout=layout, dpi=100, size=(1100, 850))
        if subtract_leading_order
            filename_extension = "_without_leading_order"
        else
            filename_extension = ""
        end
        #annotate!((0., 0., directory * filename), subplot=1)
        savefig("HSX_coil_hifi_B_vector" * filename_extension * ".pdf")
    end

end

function plot_high_fidelity_B_vector_for_HSX_coil_rectangular_xsection(
    filename="HSX_rectangular_xsection_B_a0.02_b0.01_rtol1.0e-5_atol1.0e-5_nphi4_nu25_nv26_2023-06-21T16.02.19.567.dat"
)
    directory = "/Users/mattland/Box/work23/20230621-01-rectangular_xsection_B/"
    file = open(directory * filename, "r")

    # Read and parse header:
    line = readline(file)
    line = readline(file)
    splitline = split(line, ",")
    a = parse(Float64, splitline[1])
    b = parse(Float64, splitline[2])
    I = parse(Float64, splitline[3])
    reltol = parse(Float64, splitline[4])
    abstol = parse(Float64, splitline[5])
    nϕ = parse(Int, splitline[6])
    nu = parse(Int, splitline[7])
    nv = parse(Int, splitline[8])
    println("Read nϕ=$(nϕ), nu=$(nu), nv=$(nv)")

    curve = get_curve("hsx", 1)
    #curve = CurveCircle(1.0)
    coil = CoilRectangularXSection(curve, I, a, b, FrameCentroid(curve))

    # Now read the main B data:
    high_fidelity_B = zeros(nϕ, nu, nv, 3)
    for jϕ in 1:nϕ
        for ju in 1:nu
            for jv = 1:nv
                for jxyz in 1:3
                    high_fidelity_B[jϕ, ju, jv, jxyz] = parse(Float64, readline(file))
                end
            end
        end
    end
    close(file)
    #@show high_fidelity_B

    # Set up grid of subplots
    n_plots = nϕ * 3 * 3
    n_cols = Int(ceil(0.9 * sqrt(n_plots)))
    n_rows = Int(ceil(n_plots / n_cols))
    @show n_plots, n_rows, n_cols
    xyz = "xyz"

    layout = (n_rows, n_cols)
    plots = Array{Any, 1}(undef, n_plots)
    scalefontsizes()
    scalefontsizes(0.5)

    Plots.gr_cbar_width[] = 0.005

    almost_one = 1 - (1e-6)
    u = [((ju - 1) / (nu - 1) * 2 - 1) * almost_one for ju in 1:nu]
    v = [((jv - 1) / (nv - 1) * 2 - 1) * almost_one for jv in 1:nv]
    #u2d = [u[ju] for jv in 1:nv, ju in 1:nu]
    #v2d = [v[jv] for jv in 1:nv, ju in 1:nu]
    u_a_over_2 = 0.5 * u * a
    v_b_over_2 = 0.5 * v * b

    for subtract_leading_order in [false, true]
        index = 1
        for jϕ in 1:nϕ
            ϕ = 2π * (jϕ - 1) / nϕ
            differential_arclength, curvature, torsion, position, tangent, normal, binormal = Frenet_frame(curve, ϕ)
            p, q = get_frame(coil.frame, ϕ, position, tangent, normal)
            κ1, κ2 = get_κ1_κ2(p, q, normal, curvature)

            B_regularized_plus_extra_term = B_filament_adaptive(coil, position; regularization=compute_regularization(coil))
            B_regularized_plus_extra_term += μ0 * coil.current * curvature / (8π) * (4 + 2 * log(2) + log(CoilForces.rectangular_xsection_δ(a, b))) * binormal
            
            B0 = zeros(nu, nv, 3)
            for ju in 1:nu
                for jv in 1:nv
                    B0_p, B0_q = CoilForces.rectangular_xsection_B0(coil, u[ju], v[jv])
                    B0[ju, jv, :] = B0_p * p + B0_q * q
                end
            end
        
            for jxyz in 1:3
                B_hifi = high_fidelity_B[jϕ, :, :, jxyz]

                B_analytic = zeros(nu, nv)
                for ju in 1:nu
                    for jv in 1:nv
                        Bκ_p, Bκ_q = CoilForces.rectangular_xsection_Bκ(coil, κ1, κ2, u[ju], v[jv])
                        Bκ_term = Bκ_p * p + Bκ_q * q
                        B_analytic[ju, jv] = (
                            B_regularized_plus_extra_term[jxyz] 
                            + B0[ju, jv, jxyz]
                            + Bκ_term[jxyz]
                            )
                    end
                end

                if subtract_leading_order
                    B_analytic -= B0[:, :, jxyz]
                    B_hifi -= B0[:, :, jxyz]
                end

                for variant in ["hifi", "analytic", "difference"]
                    if variant == "hifi"
                        data = B_hifi
                    elseif variant == "analytic"
                        data = B_analytic
                    elseif variant == "difference"
                        data = B_hifi - B_analytic
                    end

                    #maxB = maximum(leading_order_solution)
                    #minB = minimum(leading_order_solution)
                    #if maxB == minB
                    #    maxB += 1e-10
                    #end
                    #contour_levels = collect(range(minB, maxB, length=25))
                    #@show contour_levels

                    plots[index] = contour(u_a_over_2, v_b_over_2, data',
                        aspect_ratio = :equal,
                        #levels=contour_levels,
                    )
                    index += 1
                    title_str = @sprintf "B%s [Tesla] at ϕ=%.2f" xyz[jxyz] ϕ
                    #title!("B$(xyz[jxyz]) [Tesla] at ϕ=$(ϕ)")
                    title_str = variant * " " * title_str
                    title!(title_str)
                    xlabel!("ua/2 [meters]")
                    ylabel!("vb/2 [meters]")
                end
            end
        end
        plot(plots..., layout=layout, dpi=100, size=(2200, 1700))
        if subtract_leading_order
            filename_extension = "_without_leading_order"
        else
            filename_extension = ""
        end
        #annotate!((0., 0., directory * filename), subplot=1)
        savefig(directory * "HSX_rectangular_xsection_B" * filename_extension * ".pdf")
    end

end

function plot_high_fidelity_B_vector_for_HSX_coil_rectangular_xsection_for_paper(
    #filename="HSX_rectangular_xsection_B_a0.02_b0.01_rtol1.0e-5_atol1.0e-5_nphi4_nu25_nv26_2023-06-21T16.02.19.567.dat"
    filename="HSX_rectangular_xsection_B_a0.1296_b0.0568_rtol1.0e-5_atol1.0e-5_nphi1_nu64_nv65_2023-06-24T08.26.44.911.dat"
)
    # To generate hifi data for this function, run
    # CoilForces.save_high_fidelity_B_vector_for_HSX_coil_rectangular_xsection()
    directory = "/Users/mattland/Box/work23/20230621-01-rectangular_xsection_B/"
    file = open(directory * filename, "r")

    # Read and parse header:
    line = readline(file)
    line = readline(file)
    splitline = split(line, ",")
    a = parse(Float64, splitline[1])
    b = parse(Float64, splitline[2])
    I = parse(Float64, splitline[3])
    reltol = parse(Float64, splitline[4])
    abstol = parse(Float64, splitline[5])
    nϕ = parse(Int, splitline[6])
    nu = parse(Int, splitline[7])
    nv = parse(Int, splitline[8])
    println("Read nϕ=$(nϕ), nu=$(nu), nv=$(nv), I=$(I)")

    curve = get_curve("hsx", 1)
    #curve = CurveCircle(1.0)
    coil = CoilRectangularXSection(curve, I, a, b, FrameCentroid(curve))

    # Now read the main B data:
    high_fidelity_B = zeros(nϕ, nu, nv, 3)
    for jϕ in 1:nϕ
        for ju in 1:nu
            for jv = 1:nv
                for jxyz in 1:3
                    high_fidelity_B[jϕ, ju, jv, jxyz] = parse(Float64, readline(file))
                end
            end
        end
    end
    close(file)
    #@show high_fidelity_B

    scalefontsizes()
    scalefontsizes(0.9)

    Plots.gr_cbar_width[] = 0.015

    almost_one = 1 - (1e-6)
    u = [((ju - 1) / (nu - 1) * 2 - 1) * almost_one for ju in 1:nu]
    v = [((jv - 1) / (nv - 1) * 2 - 1) * almost_one for jv in 1:nv]
    #u2d = [u[ju] for jv in 1:nv, ju in 1:nu]
    #v2d = [v[jv] for jv in 1:nv, ju in 1:nu]
    u_a_over_2 = 0.5 * u * a
    v_b_over_2 = 0.5 * v * b

    jϕ = 1
    ϕ = 2π * (jϕ - 1) / nϕ
    differential_arclength, curvature, torsion, position, tangent, normal, binormal = Frenet_frame(curve, ϕ)
    p, q = get_frame(coil.frame, ϕ, position, tangent, normal)
    κ1, κ2 = get_κ1_κ2(p, q, normal, curvature)

    B_regularized_plus_extra_term = B_filament_adaptive(coil, position; regularization=compute_regularization(coil))
    B_regularized_plus_extra_term += μ0 * coil.current * curvature / (8π) * (4 + 2 * log(2) + log(CoilForces.rectangular_xsection_δ(a, b))) * binormal
    
    B0 = zeros(nu, nv, 3)
    for ju in 1:nu
        for jv in 1:nv
            B0_p, B0_q = CoilForces.rectangular_xsection_B0(coil, u[ju], v[jv])
            B0[ju, jv, :] = B0_p * p + B0_q * q
        end
    end

    B_hifi = high_fidelity_B[jϕ, :, :, :]
    B_analytic = zeros(nu, nv, 3)
    for jxyz in 1:3
        for ju in 1:nu
            for jv in 1:nv
                Bκ_p, Bκ_q = CoilForces.rectangular_xsection_Bκ(coil, κ1, κ2, u[ju], v[jv])
                Bκ_term = Bκ_p * p + Bκ_q * q
                B_analytic[ju, jv, jxyz] = (
                    B_regularized_plus_extra_term[jxyz] 
                    + B0[ju, jv, jxyz]
                    + Bκ_term[jxyz]
                    )
            end
        end
    end

    for variant in ["hifi", "analytic"]
        if variant == "hifi"
            data3 = B_hifi
            variant_str = "High fidelity"
        elseif variant == "analytic"
            data3 = B_analytic
            variant_str = "Reduced model"
        else
            throw("Should not get here")
        end

        for component in ["z", "abs"]
            if component == "z"
                data = data3[:, :, 3]
                data_str = "\$B_z\$"
                data_file_str = "Bz"
                contour_levels = collect(-0.6:0.05:0.6)
            elseif component == "abs"
                data = sqrt.(data3[:, :, 1].^2 + data3[:, :, 2].^2 + data3[:, :, 3].^2)
                data_str = "\$|B|\$"
                data_file_str = "modB"
                contour_levels = collect(0:0.05:0.85)
            else
                throw("Should not get here")
            end
                

            #maxB = maximum(leading_order_solution)
            #minB = minimum(leading_order_solution)
            #if maxB == minB
            #    maxB += 1e-10
            #end
            #contour_levels = collect(range(minB, maxB, length=25))
            #@show contour_levels

            #plots[index] = contour(u_a_over_2, v_b_over_2, data',
            # Factor of 15 is to scale 1e4 A to 150 kA
            #contour(v_b_over_2, u_a_over_2, data * 15,
            contour(v_b_over_2 * 100, u_a_over_2 * 100, data,
                aspect_ratio = :equal,
                levels=contour_levels,
                clims=(contour_levels[1], contour_levels[end]),
                dpi=100, 
                size=(250, 400),
                framestyle=:box
            )
            xlims!(-b/2 * 100, b/2 * 100)
            ylims!(-a/2 * 100, a/2 * 100)
            title_str = variant_str * " " * data_str * " [Tesla]"
            title!(title_str)
            ylabel!("\$ua/2\$ [cm]")
            xlabel!("\$vb/2\$ [cm]")
            savefig(directory * "20230621-01-HSX_rectangular_xsection_B_for_paper_" * variant * "_" * data_file_str * ".pdf")
        end
    end
end

function plot_high_fidelity_B_vector_for_HSX_coil_for_talk(
    filename="HSX_B_a0.01_rtol0.001_atol0.001_nphi2_nn30_nb31_2023-04-27T09.40.25.165.dat"
)
    directory = "/Users/mattland/Box/work23/20230426-01-HSX_coil_hifi_B_vector/"
    file = open(directory * filename, "r")

    # Read and parse header:
    line = readline(file)
    line = readline(file)
    splitline = split(line, ",")
    aminor = parse(Float64, splitline[1])
    I = parse(Float64, splitline[2])
    nϕ = parse(Int, splitline[5])
    nn = parse(Int, splitline[6])
    nb = parse(Int, splitline[7])
    println("Read nϕ=$(nϕ), nn=$(nn), nb=$(nb)")

    curve = get_curve("hsx", 1)
    #curve = CurveCircle(1.0)
    coil = CoilCircularXSection(curve, I, aminor)
    regularization = aminor * aminor / sqrt(ℯ)

    # Now read the main B data:
    high_fidelity_B = zeros(nϕ, nn, nb, 3)
    for jϕ in 1:nϕ
        for jn in 1:nn
            for jb = 1:nb
                for jxyz in 1:3
                    high_fidelity_B[jϕ, jn, jb, jxyz] = parse(Float64, readline(file))
                end
            end
        end
    end
    close(file)
    #@show high_fidelity_B

    # Set up grid of subplots
    #nϕ = 1
    #n_plots = nϕ * 3 * 2
    n_plots = 4 * 2
    n_cols = Int(ceil(0.9 * sqrt(n_plots)))
    n_rows = Int(ceil(n_plots / n_cols))
    @show n_plots, n_rows, n_cols
    xyz = ["Bx", "By", "Bz", "|B|"]

    layout = (n_rows, n_cols)
    plots = Array{Any, 1}(undef, n_plots)
    scalefontsizes()
    scalefontsizes(0.8)

    Plots.gr_cbar_width[] = 0.01

    uplot = [((jn - 1) / (nn - 1) * 2 - 1) * aminor for jn in 1:nn]
    vplot = [((jb - 1) / (nb - 1) * 2 - 1) * aminor for jb in 1:nb]
    u2d = [uplot[jn] for jb in 1:nb, jn in 1:nn]
    v2d = [vplot[jb] for jb in 1:nb, jn in 1:nn]
    ρ = @. sqrt(u2d^2 + v2d^2) / aminor
    θ = atan.(v2d, u2d)
    @show ρ
    @show θ

    for subtract_leading_order in [false, true]
        index = 1
        for jϕ in 1:1
            #for jϕ in 1:nϕ
                ϕ = 2π * (jϕ - 1) / nϕ
            differential_arclength, curvature, torsion, position, tangent, normal, binormal = Frenet_frame(curve, ϕ)

            B_regularized = B_filament_adaptive(coil, position; regularization=regularization)
            @show B_regularized
            
            for jxyz in 1:4
                # jxyz = 4 corresponds to |B|
                if jxyz < 4
                    leading_order_solution = μ0 * I / (2π * aminor * aminor) * (
                        -normal[jxyz] * v2d + binormal[jxyz] * u2d
                    )
                end
                for hifi in [true, false]
                    if jxyz < 4
                        if hifi
                            data = high_fidelity_B[jϕ, :, :, jxyz]'
                        else
                            data = zeros(nb, nn)
                            for jb in 1:nb
                                for jn in 1:nn
                                    data[jb, jn] = (
                                        B_regularized[jxyz] 
                                        + CoilForces.B_local(coil, curvature, normal[jxyz], binormal[jxyz], ρ[jb, jn], θ[jb, jn])
                                    )
                                end
                            end
                        end
                        if subtract_leading_order
                            data -= leading_order_solution
                        end
                    else
                        # jxyz = 4
                        if hifi
                            data = sqrt.(
                                high_fidelity_B[jϕ, :, :, 1] .* high_fidelity_B[jϕ, :, :, 1]
                                .+ high_fidelity_B[jϕ, :, :, 2] .* high_fidelity_B[jϕ, :, :, 2]
                                .+ high_fidelity_B[jϕ, :, :, 3] .* high_fidelity_B[jϕ, :, :, 3]
                                )'
                        else
                            data_squared = zeros(nb, nn)
                            for jb in 1:nb
                                for jn in 1:nn
                                    for jxyz2 = 1:3
                                        data_squared[jb, jn] += (
                                            B_regularized[jxyz2] 
                                            + CoilForces.B_local(coil, curvature, normal[jxyz2], binormal[jxyz2], ρ[jb, jn], θ[jb, jn])
                                        )^2
                                    end
                                    data = sqrt.(data_squared)
                                end
                            end
                        end
                    end

                    # Don't plot data outside the coil - make those points NaN
                    for jb in 1:nb
                        for jn in 1:nn
                            if uplot[jn]^2 + vplot[jb]^2 > aminor^2
                                data[jb, jn] = NaN
                            end
                        end
                    end

                    #maxB = maximum(leading_order_solution)
                    #minB = minimum(leading_order_solution)
                    #if maxB == minB
                    #    maxB += 1e-10
                    #end
                    #contour_levels = collect(range(minB, maxB, length=25))
                    ##@show contour_levels

                    length_factor = 100 # cm
                    plots[index] = contour(uplot * length_factor, vplot * length_factor, data,
                        aspect_ratio = :equal,
                        #levels=contour_levels,
                    )
                    index += 1
                    title_str = @sprintf "%s [Tesla] at ϕ=%.2f" xyz[jxyz] ϕ
                    #title!("B$(xyz[jxyz]) [Tesla] at ϕ=$(ϕ)")
                    if hifi
                        title_str = "HiFi " * title_str
                    else
                        title_str = "analytic " * title_str
                    end
                    title!(title_str)
                    xlabel!("x [cm]")
                    ylabel!("y [cm]")
                    nθ = 150
                    θplot = collect(range(0, 2π, length=nθ))
                    plot!(aminor * cos.(θplot) * length_factor, aminor * sin.(θplot) * length_factor, linewidth=1.5, color=:black, label=nothing)
                end
            end
        end
        plot(plots..., layout=layout, dpi=100, size=(1100, 850),
            margin=3mm,
        )
        if subtract_leading_order
            filename_extension = "_without_leading_order"
        else
            filename_extension = ""
        end
        #annotate!((0., 0., directory * filename), subplot=1)
        savefig("/Users/mattland/Box/work23/20230510-02_HSX_B" * filename_extension * ".pdf")
    end

end


function debug_stalling_B_integral(;
    aminor = 0.01,  # Minor radius of coil [meters]
    I = 1.0e4,  # Total current [Amperes]
    reltol = 1e-2,
    abstol = 1e-2,
)
    curve = get_curve("hsx", 1)
    #curve = CurveCircle(1.0)
    coil = CoilCircularXSection(curve, I, aminor)
    ϕ = π
    θ_shift=0.0
    eval_point = [1.4888650264603631, 0.24140179170087098, -0.2601603251227219]
    #B = B_finite_thickness_normalized(coil, eval_point; reltol=reltol, abstol=abstol, ϕ_shift=ϕ, θ_shift=θ_shift)

    differential_arclength, curvature, torsion, position, tangent, normal, binormal = Frenet_frame(curve, ϕ)
    @show eval_point
    @show position
    ρ_singular = 1
    singular_position = position + aminor * ρ_singular * (cos(θ_shift) * normal + sin(θ_shift) * binormal)
    @show singular_position

    maxevals = 5000
    source_points = zeros(maxevals, 3)
    B_evals = zeros(maxevals, 3)

    # Next comes a copy of B_finite_thickness_normalized()
    myindex = 0
    function Biot_savart_cubature_func(xp)
        myindex += 1
        if myindex % 10000 == 0
            println("myindex ", myindex)
        end
        source_points[myindex, :] = xp
        B_evals[myindex, :] = B_finite_thickness_integrand(coil, xp[1], xp[2], xp[3], eval_point)
        return B_evals[myindex, :]
    end

    ϕ_shift = ϕ
    Biot_savart_xmin = [0, θ_shift, ϕ_shift]
    Biot_savart_xmax = [1, θ_shift + 2π, ϕ_shift + 2π]

    val, err = hcubature(
        Biot_savart_cubature_func, 
        Biot_savart_xmin,
        Biot_savart_xmax;
        atol=abstol,
        rtol=reltol,
        maxevals=maxevals ÷ 2,
        #initdiv=10,
    )
    #print("Number of function evals: ", myindex)
    start_index = 1
    end_index = myindex
    p1 = plot(source_points[start_index:end_index, 1], ylabel="ρ", label=false)
    p2 = plot(source_points[start_index:end_index, 2], ylabel="θ", label=false)
    p3 = plot(source_points[start_index:end_index, 3], ylabel="ϕ", label=false)
    p4 = plot(abs.(B_evals[start_index:end_index, 1]), ylabel="Bx", label=false, yscale=:log10)
    p5 = plot(abs.(B_evals[start_index:end_index, 2]), ylabel="By", label=false, yscale=:log10)
    p6 = plot(abs.(B_evals[start_index:end_index, 3]), ylabel="Bz", label=false, yscale=:log10)
    #return val
    #plot(p1, p2, p3, layout=(3, 1))
    plot(p1, p2, p3, p4, p5, p6, layout=(6, 1), dpi=100, size=(700, 600))
end

function debug_stalling_force_integral()
    R0 = 1.0
    I = 1.0
    a = 0.01
    b = 0.01
    ϕ = 0.0
    reltol = 1e-3
    abstol = 1e-3
    curve = CurveCircle(R0)
    coil = CoilRectangularXSection(curve, I, a, b, FrameCircle())

    maxevals = 1000
    source_points = zeros(maxevals, 2)
    force_integrand_evals = zeros(maxevals, 3)

    # Next comes a copy of force_finite_thickness()
    dℓdϕ, κ, τ, r0, tangent, normal, binormal = Frenet_frame(coil.curve, ϕ)
    p, q = get_frame(coil.frame, ϕ, r0, tangent, normal)
    κ1, κ2 = CoilForces.get_κ1_κ2(p, q, normal, κ)
    
    myindex = 0
    function force_cubature_func(xp)
        myindex += 1
        if myindex % 200 == 0
            println("myindex ", myindex)
        end
        source_points[myindex, :] = xp
        u = xp[1]
        v = xp[2]
        u_a_over_2 = 0.5 * u * coil.a
        v_b_over_2 = 0.5 * v * coil.b
        sqrtg = 1 - u_a_over_2 * κ1 - v_b_over_2 * κ2
        r_eval = r0 + u_a_over_2 * p + v_b_over_2 * q

        B1 = B_finite_thickness_normalized(
            coil,
            r_eval,
            reltol=reltol,
            abstol=abstol,
            ϕ_shift=ϕ,
            u_range=(-1, u),
            v_range=(-1, v),
        )
        B2 = B_finite_thickness_normalized(
            coil,
            r_eval,
            reltol=reltol,
            abstol=abstol,
            ϕ_shift=ϕ,
            u_range=(u, 1),
            v_range=(-1, v),
        )
        B3 = B_finite_thickness_normalized(
            coil,
            r_eval,
            reltol=reltol,
            abstol=abstol,
            ϕ_shift=ϕ,
            u_range=(-1, u),
            v_range=(v, 1),
        )
        B4 = B_finite_thickness_normalized(
            coil,
            r_eval,
            reltol=reltol,
            abstol=abstol,
            ϕ_shift=ϕ,
            u_range=(u, 1),
            v_range=(v, 1),
        )
        to_return  = sqrtg * cross(tangent, B1 + B2 + B3 + B4)

        """
        B = B_finite_thickness_normalized(
            coil,
            r_eval,
            reltol=reltol,
            abstol=abstol,
            ϕ_shift=ϕ,
        )
        to_return  = sqrtg * cross(tangent, B)
        """

        force_integrand_evals[myindex, :] = to_return
        return to_return
    end

    force_xmin = [-1, -1]
    force_xmax = [1, 1]
    
    val, err = hcubature(
        force_cubature_func, 
        force_xmin,
        force_xmax;
        atol=abstol * 10,
        rtol=reltol * 10,
        maxevals=maxevals - 50,
    )
    #print("Number of function evals: ", myindex)

    start_index = 1
    end_index = myindex
    println("max and min of force_integrand_evals y:", maximum(force_integrand_evals[start_index:end_index, 2]), " ", minimum(force_integrand_evals[start_index:end_index, 2]))
    println("max and min of force_integrand_evals z:", maximum(force_integrand_evals[start_index:end_index, 3]), " ", minimum(force_integrand_evals[start_index:end_index, 3]))

    """
    p1 = plot(source_points[start_index:end_index, 1], ylabel="u", label=false)
    p2 = plot(source_points[start_index:end_index, 2], ylabel="v", label=false)
    p3 = plot(abs.(force_integrand_evals[start_index:end_index, 1]), ylabel="f integrand x", label=false, yscale=:log10)
    p4 = plot(abs.(force_integrand_evals[start_index:end_index, 2]), ylabel="f integrand y", label=false)
    p5 = plot(abs.(force_integrand_evals[start_index:end_index, 3]), ylabel="f integrand z", label=false)
    #return val
    #plot(p1, p2, p3, layout=(3, 1))
    plot(p1, p2, p3, p4, p5, layout=(5, 1), dpi=100, size=(700, 600))
    """
    scatter(
        source_points[start_index:end_index, 1],
        source_points[start_index:end_index, 2],
        zcolor=abs.(force_integrand_evals[start_index:end_index, 1]),
        xlabel="u",
        ylabel="v",
        dpi=100, 
        size=(700, 600),
        ms=2,
        msw=0,
        label=false,
    )
end

function reproduce_Sienas_plot_of_locally_circular_approx()
    curve = get_curve("hsx", 1)
    current = 1e6
    R_eff = curve_length(curve) / (2π)
    aminor = 0.01 * R_eff
    @show aminor
    coil = CoilCircularXSection(curve, current, aminor)
    nϕ = 200
    ϕ = [(jϕ - 1) * 2π / nϕ for jϕ in 1:nϕ]
    forces = zeros(nϕ, 3)
    for j in 1:nϕ
        forces[j, :] = force_locally_circular_approximation(coil, ϕ[j])
    end

    plot(ϕ, forces[:, 1], label=false, color=:darkorange, dpi=100, size=(510, 400))
    xlabel!("ϕ")
    title!("dFₓ/dℓ for HSX coil 1, Garren's locally circular approximation")
    ylims!(-3.7e6, 2.25e6)

end

function reproduce_Sienas_plot_of_exact_and_locally_circular_approx()
    directory = "/Users/mattland/Box/work23/20230429-01_hifi_force_for_HSX_coil_1/"
    filename = "20230429-01-hifi_force_for_HSX_coil_1_a0.0032695461797682913_nphi100_rtol0.01_atol0.01_2023-04-30T16.14.53.488.dat"
    csv_file = CSV.File(directory * filename)
    
    curve = get_curve("hsx", 1)
    current = 1e6
    R_eff = curve_length(curve) / (2π)
    aminor = 0.01 * R_eff
    @show aminor
    coil = CoilCircularXSection(curve, current, aminor)
    nϕ = 200
    ϕ = [(jϕ - 1) * 2π / nϕ for jϕ in 1:nϕ]
    forces = zeros(nϕ, 3)
    for j in 1:nϕ
        forces[j, :] = force_locally_circular_approximation(coil, ϕ[j])
    end

    # Save 1D result to a text file
    new_filename = "20230429-01-locally_circular_approx_force_for_HSX_coil_1_a$(aminor)_I$(current)_nphi$(nϕ).dat"
    open(directory * new_filename, "w") do file
        write(file, "ϕ, force_x/length, force_y/length, force_z/length\n")
        for jϕ in 1:nϕ
            write(file, "$(ϕ[jϕ]), $(forces[jϕ, 1]), $(forces[jϕ, 2]), $(forces[jϕ, 3])\n")
        end
    end

    plot(csv_file.ϕ, csv_file[" force_x/length"], label="exact", dpi=100, size=(510, 400))
    plot!(ϕ, forces[:, 1], label="locally circular approx", color=:darkorange)
    xlabel!("ϕ")
    title!("dFₓ/dℓ for HSX coil 1")
    ylims!(-3.7e6, 2.25e6)
end

function compare_hifi_force_to_1D_for_HSX()
    #filename = "20230429-01-hifi_force_for_HSX_coil_1_a0.0032695461797682913_nphi100_rtol0.01_atol0.01_2023-04-30T16.14.53.488.dat"
    #aminor = 0.0032695461797682913
    #current = 1e6

    #filename = "20230429-01-hifi_force_for_HSX_coil_1_a0.01_nphi100_rtol0.001_atol0.001_2023-04-30T07.30.37.591.dat"
    #aminor = 0.01
    #current = 1.5e5

    #filename = "20230429-01-hifi_force_for_HSX_coil_1_a0.1_nphi100_rtol0.001_atol0.001_2023-05-02T11.26.56.050.dat"
    #aminor = 0.1
    #current = 1e6

    filename = "20230429-01-hifi_force_for_HSX_coil_1_a0.001_nphi100_rtol0.0001_atol0.0001_2023-05-01T20.47.56.154.dat"
    aminor = 0.001
    current = 1e6

    directory = "/Users/mattland/Box/work23/20230429-01_hifi_force_for_HSX_coil_1/"
    csv_file = CSV.File(directory * filename)
    
    curve = get_curve("hsx", 1)
    @show aminor
    regularization = aminor * aminor / sqrt(ℯ)
    coil = CoilCircularXSection(curve, current, aminor)
    nϕ = length(csv_file.ϕ)
    ϕs = [(jϕ - 1) * 2π / nϕ for jϕ in 1:nϕ]
    forces = zeros(nϕ, 3)
    for j in 1:nϕ
        ϕ = ϕs[j]
        r_eval, tangent_vec = position_and_tangent(curve, ϕ)
        B = B_filament_adaptive(coil, r_eval; regularization=regularization)
        forces[j, :] = current * cross(tangent_vec, B)
    end

    modF = sqrt.(forces[:, 1].^2 + forces[:, 2].^2 + forces[:, 3].^2)
    @show modF
    normalization = sqrt(sum(modF.^2) / nϕ)
    @show normalization

    p1 = plot(csv_file.ϕ, csv_file[" force_x/length"], label="high fidelity")
    plot!(ϕs, forces[:, 1], linestyle=:dash, label="1D filament model", color=:darkorange)
    xlabel!("ϕ")
    title!("dFx/dℓ")

    p2 = plot(csv_file.ϕ, csv_file[" force_y/length"], label="high fidelity")
    plot!(ϕs, forces[:, 2], linestyle=:dash, label="1D filament model", color=:darkorange)
    xlabel!("ϕ")
    title!("dFy/dℓ")

    p3 = plot(csv_file.ϕ, csv_file[" force_z/length"], label="high fidelity")
    plot!(ϕs, forces[:, 3], linestyle=:dash, label="1D filament model", color=:darkorange)
    xlabel!("ϕ")
    title!("dFz/dℓ")

    p4 = plot(ϕs, (csv_file[" force_x/length"] - forces[:, 1]) / normalization, label=false)
    xlabel!("ϕ")
    title!("Relative difference in dFx/dℓ")

    p5 = plot(ϕs, (csv_file[" force_y/length"] - forces[:, 2]) / normalization, label=false)
    xlabel!("ϕ")
    title!("Relative difference in dFy/dℓ")

    p6 = plot(ϕs, (csv_file[" force_z/length"] - forces[:, 3]) / normalization, label=false)
    xlabel!("ϕ")
    title!("Relative difference in dFz/dℓ")

    plot(p1, p2, p3, p4, p5, p6, layout=(2, 3), dpi=100, size=(1100, 850), 
        plot_title="HSX coil 1\n$(filename)",
    )
    directory2 = "/Users/mattland/Box/work23/"
    savefig(directory2 * "HSX_coil_force_hifi_vs_filament_a$(aminor).pdf")

    # Save 1D result to a text file
    new_filename = "20230429-01-regularized_1D_BS_force_for_HSX_coil_1_a$(aminor)_I$(current)_nphi$(nϕ).dat"
    open(directory * new_filename, "w") do file
        write(file, "ϕ, force_x/length, force_y/length, force_z/length\n")
        for jϕ in 1:nϕ
            write(file, "$(ϕs[jϕ]), $(forces[jϕ, 1]), $(forces[jϕ, 2]), $(forces[jϕ, 3])\n")
        end
    end
    println("Finished.")
end


function plot_hifi_force_and_1D_for_HSX_for_talk()
    #filename = "20230429-01-hifi_force_for_HSX_coil_1_a0.0032695461797682913_nphi100_rtol0.01_atol0.01_2023-04-30T16.14.53.488.dat"
    #aminor = 0.0032695461797682913
    #current = 1e6

    #filename = "20230429-01-hifi_force_for_HSX_coil_1_a0.01_nphi100_rtol0.001_atol0.001_2023-04-30T07.30.37.591.dat"
    #aminor = 0.01
    #current = 1.5e5

    filenames = [
        "20230429-01-hifi_force_for_HSX_coil_1_a0.1_nphi100_rtol0.001_atol0.001_2023-05-02T11.26.56.050.dat",
        "20230429-01-hifi_force_for_HSX_coil_1_a0.01_nphi100_rtol0.001_atol0.001_2023-04-30T07.30.37.591.dat",
        #"20230429-01-hifi_force_for_HSX_coil_1_a0.0032695461797682913_nphi100_rtol0.01_atol0.01_2023-04-30T16.14.53.488.dat",
    ]
    aminors = [
        0.1,
        0.01,
        #0.0032695461797682913,
    ]
    currents = [
        1e6,
        150e3,
        1e6,
    ]

    # Rescale results to a current of 150 kA
    current_scale_factors = (150e3) ./ currents
    #current_scale_factors = [1, 1]

    #filename = "20230429-01-hifi_force_for_HSX_coil_1_a0.001_nphi100_rtol0.0001_atol0.0001_2023-05-01T20.47.56.154.dat"
    #aminor = 0.001
    #current = 1e6

    directory = "/Users/mattland/Box/work23/20230429-01_hifi_force_for_HSX_coil_1/"

    plots = Vector{Any}(undef, 2)
    scalefontsizes()
    #scalefontsizes(0.5)

    for ja in 1:2
        filename = filenames[ja]
        aminor = aminors[ja]
        current = currents[ja]
        csv_file = CSV.File(directory * filename)
        
        curve = get_curve("hsx", 1)
        @show aminor
        regularization = aminor * aminor / sqrt(ℯ)
        coil = CoilCircularXSection(curve, current, aminor)
        nϕ = length(csv_file.ϕ)
        ϕs = [(jϕ - 1) * 2π / nϕ for jϕ in 1:nϕ]
        forces = zeros(nϕ, 3)
        for j in 1:nϕ
            ϕ = ϕs[j]
            r_eval, tangent_vec = position_and_tangent(curve, ϕ)
            B = B_filament_adaptive(coil, r_eval; regularization=regularization)
            forces[j, :] = current * cross(tangent_vec, B)
        end

        # Note that force ∝ current^2:
        normalization = 1e3 / (current_scale_factors[ja]^2)
        plots[ja] = plot(
            csv_file.ϕ, csv_file[" force_x/length"] / normalization, label=false, framestyle = :box,
            linewidth=2,
        )
        plot!(
            ϕs, forces[:, 1] / normalization, linestyle=:dash, label=false, color=:darkorange,
            linewidth=3
        )
        xlabel!("ϕ")
        title!("a = $(Int(aminor * 100)) cm")
        xlims!(0, 2π)
    end

    plot(plots..., layout=(1, 2), dpi=100, size=(950, 400), 
        #plot_title="HSX coil 1",
    )
    savefig("/Users/mattland/Box/work23/20230508-01-HSX_coil_force_hifi_vs_filament.pdf")

    println("Finished.")
end

function save_inductance_a_scan()
    reltol = 1e-6
    abstol = 1e-6

    hsx = true
    #hsx = false

    aminors = 10 .^ collect(((-4.0):(0.0625):(-1.0625)))
    #aminors = 10 .^ collect(((-3.0):(0.0625):(0)))
    #aminors = 10 .^ collect(((-2.0):(0.5):(0)))
    println("Values of a/R that will be evaluated: ", aminors)

    if hsx
        curve = get_curve("hsx", 1)
        config_str = "hsx"
    else
        # Major radius of coil [meters]
        R0 = 1.0
        curve = CurveCircle(R0)
        config_str = "circle"
    end

    # Total current [Amperes]
    I = 1.0

    inductances_hifi = similar(aminors)
    inductances_filament = similar(aminors)
    times = similar(aminors)
    for ja in eachindex(aminors)
        a = aminors[ja]
        println("a = ", a)
        coil = CoilCircularXSection(curve, I, a)

        @time L_filament = inductance_filament_adaptive(coil; abstol=0, reltol=1e-9)
        time_data = @timed L_hifi = inductance_finite_thickness(coil; reltol=reltol, abstol=abstol)
        times[ja] = time_data.time
        inductances_hifi[ja] = L_hifi
        inductances_filament[ja] = L_filament
        println("  time: $(time_data.time)  L_filament: $(L_filament)  L_hifi: $(L_hifi)  ratio: $(L_filament / L_hifi)")
    end

    directory = "/Users/mattland/Box/work23/20230517-01-inductance_a_scans/"
    date_str = replace("$(Dates.now())", ":" => ".")
    filename = "inductance_$(config_str)_rtol_$(reltol)_atol_$(abstol)_$(date_str).dat"
    open(directory * filename, "w") do file
        write(file, "a,L_hifi,L_filament,time\n")
        for ja in eachindex(aminors)
            write(file, "$(aminors[ja]),$(inductances_hifi[ja]),$(inductances_filament[ja]),$(times[ja])\n")
        end
    end
    
end

function plot_inductance_a_scan()
    directory = "/Users/mattland/Box/work23/20230517-01-inductance_a_scans/"
    filenames = [
        "inductance_circle_rtol_0.001_atol_0.001_2023-05-17T18.07.45.692.dat",
        "inductance_circle_rtol_0.0001_atol_0.0001_2023-05-17T18.08.57.739.dat",
        "inductance_circle_rtol_1.0e-5_atol_1.0e-5_2023-05-17T18.16.35.204.dat",
        "inductance_circle_rtol_1.0e-6_atol_1.0e-6_2023-05-17T19.04.53.314.dat",
    ]
    filenames = [
        "inductance_hsx_rtol_0.01_atol_0.01_2023-05-17T21.11.27.941.dat",
        "inductance_hsx_rtol_0.001_atol_0.001_2023-05-18T02.29.30.339.dat",
    ]
    filenames = [
        "inductance_hsx_rtol_0.01_atol_0.01_2023-05-19T06.06.42.412.dat",
        "inductance_hsx_rtol_0.001_atol_0.001_2023-05-19T06.14.05.159.dat",
        "inductance_hsx_rtol_0.0001_atol_0.0001_2023-05-19T06.36.41.152.dat",
        "inductance_hsx_rtol_1.0e-5_atol_1.0e-5_2023-05-23T13.47.07.432.dat",
        "inductance_hsx_rtol_1.0e-6_atol_1.0e-6_2023-05-27T08.08.26.215.dat",
    ]
    n = length(filenames)
    plot()
    #    xscale=:log10,
    #    yscale=:log10,
    #)
    for j in 1:n
        filename = filenames[j]
        f = CSV.File(directory * filename)
        data = @. (f.L_hifi - f.L_filament) / f.L_hifi
    
        plot!(f.a, abs.(data), 
            xscale=:log10,
            yscale=:log10,
            label=filename,
            minorgrid=true,
        )
    end
    xlabel!("aminor")
    title!("|L_hifi - L_filament| / L_hifi")
end

function plot_inductance_convergence()
    curve = get_curve("hsx", 1)
    regularization = 0.01 ^ 2
    n = 300

    # Generate numbers of quadrature points to try:
    nns = 30
    ns = [Int(round(10 ^ x)) for x in range(1.0, 3.0, length=nns)]
    
    inductance_results = zeros(nns)
    inductance_results_singularity_subtraction = zeros(nns)
    for jn in 1:nns
        n = ns[jn]
        inductance_results[jn] = CoilForces.inductance_filament_fixed(curve, regularization, n)
        inductance_results_singularity_subtraction[jn] = CoilForces.inductance_filament_fixed_singularity_subtraction(curve, regularization, n)
    end
    @show inductance_results
    @show inductance_results_singularity_subtraction

    scatter(
        ns,
        inductance_results,
        xscale=:log10,
        #yscale=:log10,
        minorgrid=true,
        label="original"
    )
    scatter!(
        ns,
        inductance_results_singularity_subtraction,
        label="singularity subtraction"
    )
end

function save_inductance_b_scan_rectangular_xsection(;
    a=0.01,
    hsx=true,
    reltol = 1e-1,
    abstol = 1e-1,
)

    #bs = a * 10 .^ collect(((-1.0):(0.125):(1)))
    bs = a * 10 .^ collect(((-1.0):(0.125):(1)))
    nb = length(bs)
    println("a = $(a)")
    println("Values of b that will be evaluated: ", bs)

    # Total current [Amperes]
    I = 1.0

    if hsx
        config_str = "hsx"
    else
        config_str = "circle"
    end

    inductances_hifi = similar(bs)
    inductances_filament = similar(bs)
    times = similar(bs)
    println("Number of threads: ", Threads.nthreads())
    Threads.@threads for jb in 1:nb
        # CurveXYZFourier has buffers that need to have distinct contents in
        # different threads, so we define the curve here inside the loop, where
        # all new variables are distinct for each thread.
        if hsx
            curve = get_curve("hsx", 1)
        else
            # Major radius of coil [meters]
            R0 = 1.0
            curve = CurveCircle(R0)
        end
    
        b = bs[jb]
        #println("b = ", b)
        println("Thread $(Threads.threadid()) is handling jb = $(jb) of $(nb): b = ", b)
        coil = CoilRectangularXSection(curve, I, a, b, FrameCentroid(curve))

        @time L_filament = inductance_filament_adaptive(coil; abstol=0, reltol=1e-3)
        #if sqrt(a^2 + b^2) < 1 / 11.5
        #if sqrt(a^2 + b^2) < 0.075
        if true
            time_data = @timed L_hifi = inductance_finite_thickness(coil; reltol=reltol, abstol=abstol)
        else
            time_data = @timed L_hifi = NaN
        end
        times[jb] = time_data.time
        inductances_hifi[jb] = L_hifi
        inductances_filament[jb] = L_filament
        println("  time: $(time_data.time)  L_filament: $(L_filament)  L_hifi: $(L_hifi)  ratio: $(L_filament / L_hifi)")
    end

    directory = "/Users/mattland/Box/work23/20230517-01-inductance_a_scans/"
    date_str = replace("$(Dates.now())", ":" => ".")
    filename = "inductance_rectangular_xsection_$(config_str)_a$(a)_rtol$(reltol)_atol$(abstol)_$(date_str).dat"
    open(directory * filename, "w") do file
        write(file, "a, reltol, abstol\n")
        write(file, "$(a),$(reltol),$(abstol)\n")
        write(file, "b,L_hifi,L_filament,time\n")
        for jb in eachindex(bs)
            write(file, "$(bs[jb]),$(inductances_hifi[jb]),$(inductances_filament[jb]),$(times[jb])\n")
        end
    end
    
end

function plot_inductance_b_scan_rectangular_xsection()
    directory = "/Users/mattland/Box/work23/20230517-01-inductance_a_scans/"
    filenames = [
        "inductance_rectangular_xsection_hsx_a0.001_rtol0.1_atol0.1_2023-06-22T12.55.18.595.dat",
        "inductance_rectangular_xsection_hsx_a0.01_rtol0.1_atol0.1_2023-06-22T16.53.30.360.dat",
    ]
    filenames = [
        "inductance_rectangular_xsection_hsx_a0.001_rtol0.001_atol0.001_2023-06-27T06.22.45.637.dat",
        "inductance_rectangular_xsection_hsx_a0.001_rtol0.0001_atol0.0001_2023-06-27T06.28.31.044.dat",
        "inductance_rectangular_xsection_hsx_a0.01_rtol0.001_atol0.001_2023-06-26T19.22.27.950.dat",
        "inductance_rectangular_xsection_hsx_a0.01_rtol0.0001_atol0.0001_2023-06-26T19.31.00.200.dat",
        "inductance_rectangular_xsection_hsx_a0.1_rtol0.001_atol0.001_2023-06-27T05.30.01.425.dat",
        "inductance_rectangular_xsection_hsx_a0.1_rtol0.0001_atol0.0001_2023-06-27T06.43.02.433.dat",
    ]
    filenames = [
        "inductance_rectangular_xsection_hsx_a0.1_rtol0.001_atol0.001_2023-06-27T05.30.01.425.dat",
        "inductance_rectangular_xsection_hsx_a0.1_rtol0.0001_atol0.0001_2023-06-27T06.43.02.433.dat",
    ]
    filenames = [
        "inductance_rectangular_xsection_hsx_a0.001_rtol0.0001_atol0.0001_2023-06-27T06.28.31.044.dat",
        "inductance_rectangular_xsection_hsx_a0.01_rtol0.0001_atol0.0001_2023-06-26T19.31.00.200.dat",
        "inductance_rectangular_xsection_hsx_a0.1_rtol0.0001_atol0.0001_2023-06-27T06.43.02.433.dat",
    ]
    n = length(filenames)

    scalefontsizes()
    scalefontsizes(1.0)
    plot(size=(500, 500))
    #    xscale=:log10,
    #    yscale=:log10,
    #)
    # For colors, see http://juliagraphics.github.io/Colors.jl/stable/namedcolors/
    colors_filament = [:salmon, :limegreen, :deepskyblue]  # light
    colors_hifi = [:firebrick, :darkgreen, :blue3] # dark
    for j in 1:n
        filename = filenames[j]
        f = CSV.File(directory * filename, header=3)
    
        plot!(f.b, f.L_hifi * 1e6, 
            xscale=:log10,
            minorgrid=true,
            label=false,
            lw=4,
            color=colors_hifi[j]
        )
        plot!(f.b, 
            f.L_filament * 1e6,
            label=false,
            lw=2,
            ls=:dot,
            color=colors_filament[j]
        )
    end
    xlabel!("Conductor thickness b [meters]")
    ylabel!("Self-inductance [\$\\mu\$H]")
    ylims!(0, 3)
    """
    annotate!(5e-4, 2.4, text("a=0.001 m,\n2D", colors_filament[1]))
    annotate!(6e-3, 1.4, text("a=0.01 m,\n2D", colors_filament[2]))
    annotate!(0.15, 0.3, text("a=0.1 m,\n2D", colors_filament[3]))

    annotate!(4e-3, 2.7, text("a=0.001 m,\n6D", colors_hifi[1]))
    annotate!(8e-2, 1.6, text("a=0.01 m,\n6D", colors_hifi[2]))
    annotate!(4.5e-1, 0.75, text("a=0.1 m,\n6D", colors_hifi[3]))
    """
    filament_string = "2D"
    annotate!(5e-4, 2.4, text("a=0.001 m,\n" * filament_string, colors_filament[1]))
    annotate!(6e-3, 1.4, text("a=0.01 m,\n" * filament_string, colors_filament[2]))
    annotate!(0.15, 0.3, text("a=0.1 m,\n" * filament_string, colors_filament[3]))

    hifi_string = "6D"
    annotate!(4e-3, 2.7, text("a=0.001 m,\n" * hifi_string, colors_hifi[1]))
    annotate!(8e-2, 1.6, text("a=0.01 m,\n" * hifi_string, colors_hifi[2]))
    annotate!(4.5e-1, 0.75, text("a=0.1 m,\n" * hifi_string, colors_hifi[3]))

    savefig(directory * "20230517-01-inductance_rectangular_xsection_hsx_ab_scan.pdf")
end

function save_HSX_rectangular_coil_shape()
    curve = get_curve("hsx", 1)

    # Dimensions in Singh paper:
    a = 0.12
    b = 0.06

    # Dimensions in David Anderson's note:
    a = 0.1296
    b = 0.0568

    coil = CoilRectangularXSection(curve, 1.0, a, b, FrameCentroid(curve))
    regularization = compute_regularization(coil)
    n = 200  # Number of points
    filename = "/Users/mattland/Box/work23/20230614-01-hsx_rectangular_coil"
    CoilForces.save(coil, filename * ".dat", n)

    # Now save a second file with the self-force
    open(filename * "_force.dat", "w") do file
        write(file, "Fx, Fy, Fz\n")
        for j in 1:n
            ϕ = 2π * (j - 1) / n
            differential_arclength, curvature, torsion, position, tangent, normal, binormal = Frenet_frame(curve, ϕ)
            B = B_filament_adaptive(coil, position; regularization=regularization)
            F = cross(tangent, B)
            write(file, "$(F[1]), $(F[2]), $(F[3])\n")
        end
    end
end