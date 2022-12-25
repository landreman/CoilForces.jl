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

function plot_force_convergence()
    coil_num = 2

    # point at which to evaluate the force:
    ϕ0 = 0.0

    curve = get_curve("hsx", coil_num)
    # curve = CurveCircle(2.2)

    current = -1.5e5

    # minor radius of conductor:
    a = 0.01

    coil = Coil(curve, current, a)
    δ = a * a / sqrt(ℯ)

    # Generate numbers of quadrature points to try:
    nns = 20
    ns = [Int(round(10 ^ x)) for x in range(2.0, 4.0, length=nns)]

    r_eval = γ(curve, ϕ0)
    force_per_unit_length = zeros(nns)
    for jn in 1:nns
        B = B_filament_fixed(coil, r_eval, ns[jn], regularization=δ)
        force_per_unit_length[jn] = current * B[3]
    end

    scatter(ns, force_per_unit_length, xscale=:log10)
    xlabel!("number of quadrature points")
    ylabel!("Force per unit length [N / m]")
end