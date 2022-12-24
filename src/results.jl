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

    using Plots
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
    ϕ0 = -2.0

    coil = Coil(curve, current, a)
    δ = a * a / sqrt(ℯ)

    nϕ = 1000

    #ϕp = (collect(1:nϕ) .- 1) * 2π / nϕ .- π
    ϕp = collect(range(-2.01, -1.99, length=nϕ))
    #ϕp = collect(range(-π, π, length=nϕ))
    integrand = zeros(nϕ, 3)
    smoothed_integrand = zeros(nϕ, 3)
    r_eval = γ(curve, ϕ0)
    r_prime = dγdϕ(curve, ϕ0)
    r_prime_prime = d2γdϕ2(curve, ϕ0)
    dℓdϕ_squared = normsq(r_prime)
    for j in 1:nϕ
        integrand[j, :] = d_B_d_ϕ(coil, ϕp[j], r_eval, regularization=δ)
        smoothed_integrand[j, :] = (integrand[j, :] + μ0 * current / (4π) * 0.5 * cross(r_prime_prime, r_prime)
            * (2 - 2 * cos(ϕp[j] - ϕ0)) / (((2 - 2 * cos(ϕp[j] - ϕ0)) * dℓdϕ_squared + δ) ^ 1.5))
    end

    using Plots
    plot()
    xyz = "xyz"
    #plot(ϕ, integrand[:, 1], label="x")
    #plot!(ϕ, integrand[:, 2], label="y")
    #plot!(ϕ, integrand[:, 3], label="z")
    for j in 1:3
        #plot!(ϕp, integrand[:, j], label=xyz[j:j])
        plot!(ϕp, smoothed_integrand[:, j], label="smoothed, " * xyz[j:j])
    end
    xlabel!("Coil parameter ϕ")
    title!("Regularized Biot-Savart integrand for HSX coil $(coil_num) at ϕ=$(ϕ0)")

end
