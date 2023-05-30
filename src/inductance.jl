function analytic_inductance_for_circular_coil(coil::CoilCircularXSection)
    # Assert that curve type is a CurveCircle:
    coil.curve::CurveCircle
    R = coil.curve.R0
    a = coil.aminor
    return μ0 * R * (log(8 * R / a) - 1.75)
end

"""
    inductance_filament_integrand(coil::CoilCircularXSection, regularization, ϕ, ϕp)

Integrand for calculating the self-inductance using the regularized filament
method. Note that the prefactor of μ0 / (4π) is not included.
"""
function inductance_filament_integrand(curve, regularization, ϕ, ϕp)
    data1 = γ_and_derivative(curve, ϕ)
    data2 = γ_and_derivative(curve, ϕp)

    position1 = @view data1[:, 1]
    position2 = @view data2[:, 1]
    r_prime1 = @view data1[:, 2]
    r_prime2 = @view data2[:, 2]

    return (
        dot(r_prime1, r_prime2) / sqrt(
            normsq(position1 - position2) + regularization
        )
    )
end

"""
    inductance_filament_adaptive(coil::CoilCircularXSection; reltol=1e-8, abstol=1e-14)

Evaluate the Biot-Savart law for a coil in the approximation that the coil is an
infinitesmally thin filament. Use adaptive quadrature.
"""
function inductance_filament_adaptive(coil::CoilCircularXSection; reltol=1e-8, abstol=1e-14)
    aminor = coil.aminor
    regularization = aminor * aminor / √ℯ

    function inductance_filament_integrand_wrapper(ϕs) 
        # Shift 2nd angle so the near-singularity is at the boundary of the
        # integration domain. This improves speed significantly.
        return inductance_filament_integrand(coil.curve, regularization, ϕs[1], ϕs[2] + ϕs[1])
    end

    val, err = hcubature(
        inductance_filament_integrand_wrapper,
        [0, 0],
        [2π, 2π];
        rtol=reltol,
        atol=abstol,
    )
    return val * μ0 / (4π)
end

"""
Compute the self-inductance of a coil via a 6D integral, accounting for the
finite thickness.
"""
function inductance_finite_thickness(coil::CoilCircularXSection; reltol=1e-3, abstol=1e-5)

    function inductance_cubature_func(xp)
        ρ = xp[1]
        θ = xp[2]
        ϕ = xp[3]
        dℓdϕ, κ, τ, r0, tangent, normal, binormal = Frenet_frame(coil.curve, ϕ)
        r = ρ * coil.aminor
        sinθ, cosθ = sincos(θ)
        sqrtg = (1 - κ * r * cosθ) * ρ * dℓdϕ
        r_eval = r0 + r * cosθ * normal + r * sinθ * binormal
        A = A_finite_thickness_normalized(
            coil,
            r_eval,
            reltol=reltol,
            abstol=abstol,
            ϕ_shift=ϕ,
            θ_shift=θ,
        )
        return sqrtg * dot(tangent, A)
    end

    inductance_xmin = [0, 0, 0]
    inductance_xmax = [1, 2π, 2π]
    
    val, err = hcubature(
        inductance_cubature_func, 
        inductance_xmin,
        inductance_xmax;
        atol=abstol,
        rtol=reltol
    )
    A_prefactor = μ0 * coil.current / (4 * π * π)
    L_prefactor = 1 / (π * coil.current)
    return A_prefactor * L_prefactor * val
end