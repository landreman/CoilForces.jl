function analytic_inductance_for_circular_coil(coil::Coil)
    # Assert that curve type is a CurveCircle:
    coil.curve::CurveCircle
    R = coil.curve.R0
    a = coil.aminor
    return μ0 * R * (log(8 * R / a) - 1.75)
end

"""
    inductance_filament_integrand(coil::Coil, regularization, ϕ, ϕp)

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
    inductance_filament_adaptive(coil::Coil; reltol=1e-8, abstol=1e-14)

Evaluate the Biot-Savart law for a coil in the approximation that the coil is an
infinitesmally thin filament. Use adaptive quadrature.
"""
function inductance_filament_adaptive(coil::Coil; reltol=1e-8, abstol=1e-14)
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