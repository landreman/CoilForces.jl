function analytic_inductance_for_circular_coil(coil::CoilCircularXSection)
    # Assert that curve type is a CurveCircle:
    coil.curve::CurveCircle
    R = coil.curve.R0
    a = coil.aminor
    return μ0 * R * (log(8 * R / a) - 1.75)
end

function analytic_inductance_for_circular_coil(coil::CoilRectangularXSection)
    # Assert that curve type is a CurveCircle:
    coil.curve::CurveCircle
    R = coil.curve.R0
    a = coil.a
    b = coil.b
    
    # For the formula that follows, see
    # 20230531-01 Self-inductance for coils with rectangular cross-section.lyx
    return (
        μ0 * R * (log(8 * R / sqrt(a * b)) 
            + (1.0 / 12)
            + (a^4 - 6 * a^2 * b^2 + b^4) / (12 * a^2 * b^2) * log(a / b + b / a)
            - b * b / (12 * a * a) * log(b / a)
            - a * a / (12 * b * b) * log(a / b)
            - (2.0 * b) / (3 * a) * atan(a / b)
            - (2.0 * a) / (3 * b) * atan(b / a)
    ))
    
end

"""
    inductance_filament_integrand(coil::CoilCircularXSection, regularization, ϕ, ϕp)

Integrand for calculating the self-inductance using the regularized filament
method. Note that the prefactor of μ0 / (4π) is not included.
"""
function inductance_filament_integrand(curve::Curve, regularization, ϕ, ϕp)
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
    inductance_filament_integrand_singularity_subtraction(curve::Curve, regularization, ϕ, ϕp)

Integrand for calculating the self-inductance using the regularized filament
method. Note that the prefactor of μ0 / (4π) is not included.
"""
function inductance_filament_integrand_singularity_subtraction(curve::Curve, regularization, ϕ, ϕp)
    data1 = γ_and_derivative(curve, ϕ)
    data2 = γ_and_derivative(curve, ϕp)

    position1 = @view data1[:, 1]
    position2 = @view data2[:, 1]
    r_prime1 = @view data1[:, 2]
    r_prime2 = @view data2[:, 2]
    norm_sq_1 = normsq(r_prime1)

    return (
        dot(r_prime1, r_prime2) / sqrt(
            normsq(position1 - position2) + regularization
        )
        - norm_sq_1 / sqrt(
            (2 - 2 * cos(ϕ - ϕp)) * norm_sq_1 + regularization
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
    inductance_filament_fixed(coil, n, regularization)

Evaluate the self-inductance for a coil in the approximation that the coil is an
infinitesmally thin filament. Use a fixed number of quadrature points. The
pre-factor μ0 / (4π) is included.
"""
function inductance_filament_fixed(curve, regularization, n)
    dϕ = 2π / n
    val = 0.0
    for j1 in 1:n
        ϕ = (j1 - 1) * 2π / n
        for j2 in 1:n
            ϕp = (j2 - 1) * 2π / n
            val += inductance_filament_integrand(curve, regularization, ϕ, ϕp)
        end
    end
    return val * μ0 / (4π) * dϕ * dϕ
end

"""
    inductance_filament_fixed(coil, n, regularization)

Evaluate the self-inductance for a coil in the approximation that the coil is an
infinitesmally thin filament. Use a fixed number of quadrature points. The trick
of subtracting and adding a term with similar singularity is used. The
pre-factor μ0 / (4π) is included.
"""
function inductance_filament_fixed_singularity_subtraction(curve, regularization, n)
    dϕ = 2π / n

    # First, evaluate the double integral:
    val2 = 0.0
    for j1 in 1:n
        ϕ = (j1 - 1) * 2π / n
        for j2 in 1:n
            ϕp = (j2 - 1) * 2π / n
            val2 += inductance_filament_integrand_singularity_subtraction(curve, regularization, ϕ, ϕp)
        end
    end

    # Now evaluate the single integral:
    val1 = 0.0
    for j in 1:n
        ϕ = (j - 1) * 2π / n
        data = γ_and_derivative(curve, ϕ)
        r_prime = @view data[:, 2]
        r_prime_norm = norm(r_prime)
        val1 += r_prime_norm * log(64 * r_prime_norm * r_prime_norm / regularization)
    end

    return val1 * μ0 / (4π) * dϕ + val2 * μ0 / (4π) * dϕ * dϕ
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

"""
Compute the self-inductance of a coil via a 6D integral, accounting for the
finite thickness.
"""
function inductance_finite_thickness(coil::CoilRectangularXSection; reltol=1e-3, abstol=1e-5)

    function inductance_cubature_func(xp)
        u = xp[1]
        v = xp[2]
        ϕ = xp[3]
        u_a_over_2 = 0.5 * u * coil.a
        v_b_over_2 = 0.5 * v * coil.b
        α = get_winding_pack_angle(coil.winding_pack_angle, ϕ)
        sinα, cosα = sincos(α)
        dℓdϕ, κ, τ, r0, tangent, normal, binormal = Frenet_frame(coil.curve, ϕ)
        sqrtg = (1 + κ * (v_b_over_2 * sinα - u_a_over_2 * cosα)) * dℓdϕ
        r_eval = (r0
            + (u_a_over_2 * cosα - v_b_over_2 * sinα) * normal 
            + (u_a_over_2 * sinα + v_b_over_2 * cosα) * binormal
        )
        A = A_finite_thickness_normalized(
            coil,
            r_eval,
            reltol=reltol,
            abstol=abstol,
            ϕ_shift=ϕ,
        )
        return sqrtg * dot(tangent, A)
    end

    inductance_xmin = [-1, -1, 0]
    inductance_xmax = [1, 1, 2π]
    
    val, err = hcubature(
        inductance_cubature_func, 
        inductance_xmin,
        inductance_xmax;
        atol=abstol,
        rtol=reltol
    )
    A_prefactor = μ0 * coil.current / (16 * π)
    L_prefactor = 1 / (4 * coil.current)
    return A_prefactor * L_prefactor * val
end