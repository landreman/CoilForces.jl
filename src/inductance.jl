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
            - CoilForces.rectangular_xsection_k(a, b) / 2
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
function inductance_filament_adaptive(coil::Coil; reltol=1e-8, abstol=1e-14)
    regularization = compute_regularization(coil)

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
Compute the self-inductance of a coil via a 3D integral of the vector potential, accounting for the
finite thickness. This version of the function is less reliable than
inductance_finite_thickness(), but is kept here for testing.
"""
function inductance_finite_thickness_from_A(coil::CoilRectangularXSection; reltol=1e-3, abstol=1e-5)
    A_reltol = reltol * 0.1
    A_abstol = abstol * 0.1

    function inductance_cubature_func(xp)
        u = xp[1]
        v = xp[2]
        ϕ = xp[3]
        u_a_over_2 = 0.5 * u * coil.a
        v_b_over_2 = 0.5 * v * coil.b
        dℓdϕ, κ, r0, tangent, normal = Frenet_frame_without_torsion(coil.curve, ϕ)
        p, q = get_frame(coil.frame, ϕ, r0, tangent, normal)
        κ1, κ2 = CoilForces.get_κ1_κ2(p, q, normal, κ)
        #sqrtg = (0.25 * coil.a * coil.b * dℓdϕ 
        sqrtg = (dℓdϕ 
            * (1 - u_a_over_2 * κ1 - v_b_over_2 * κ2))
        r_eval = @. r0 + u_a_over_2 * p + v_b_over_2 * q
        A = A_finite_thickness_normalized(
            coil,
            r_eval,
            reltol=A_reltol,
            abstol=A_abstol,
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

"""
Compute the self-inductance of a coil via a 6D integral, accounting for the
finite thickness. This version of the function is more reliable than
"""
function inductance_finite_thickness(coil::CoilRectangularXSection; reltol=1e-3, abstol=1e-5)

    function inductance_cubature_func_cases(xp, u_case, v_case)
        u = xp[1]
        v = xp[2]
        ϕ = xp[3]
        uhat = xp[4]
        vhat = xp[5]
        ϕp = xp[6] + ϕ
        # For the change of variables that follows, see
        # 20230626-02 Moving singularity to a corner.lyx
        if u_case
            up = u + (1 - u) * uhat
            d_utilde_d_uhat = 1 - u
        else
            up = -1 + (u + 1) * uhat
            d_utilde_d_uhat = u + 1
        end
        if v_case
            vp = v + (1 - v) * vhat
            d_vtilde_d_vhat = 1 - v
        else
            vp = -1 + (v + 1) * vhat
            d_vtilde_d_vhat = v + 1
        end

        u_a_over_2 = 0.5 * u * coil.a
        v_b_over_2 = 0.5 * v * coil.b
        dℓdϕ, κ, rc, tangent, normal = Frenet_frame_without_torsion(coil.curve, ϕ)
        p, q = get_frame(coil.frame, ϕ, rc, tangent, normal)
        κ1, κ2 = CoilForces.get_κ1_κ2(p, q, normal, κ)

        up_a_over_2 = 0.5 * up * coil.a
        vp_b_over_2 = 0.5 * vp * coil.b
        dℓdϕp, κp, rcp, tangentp, normalp = Frenet_frame_without_torsion(coil.curve, ϕp)
        pp, qp = get_frame(coil.frame, ϕp, rcp, tangentp, normalp)
        κ1p, κ2p = CoilForces.get_κ1_κ2(pp, qp, normalp, κp)
        
        # Note that factors of a * b / 4 are removed from each sqrt(g) here.
        # These factors are accounted for later in "prefactor"
        sqrtg = dℓdϕ * (1 - u_a_over_2 * κ1 - v_b_over_2 * κ2)
        sqrtgp = dℓdϕp * (1 - up_a_over_2 * κ1p - vp_b_over_2 * κ2p)
        
        # This next line stores the vector r - r' in "rc"
        @. rc += u_a_over_2 * p + v_b_over_2 * q - rcp - up_a_over_2 * pp - vp_b_over_2 * qp
        return d_utilde_d_uhat * d_vtilde_d_vhat * sqrtg * sqrtgp * dot(tangent, tangentp) / (norm(rc) + 1.0e-30)
    end

    inductance_cubature_func_1(xp) = inductance_cubature_func_cases(xp, false, false)
    inductance_cubature_func_2(xp) = inductance_cubature_func_cases(xp, false, true)
    inductance_cubature_func_3(xp) = inductance_cubature_func_cases(xp, true, false)
    inductance_cubature_func_4(xp) = inductance_cubature_func_cases(xp, true, true)

    inductance_xmin = [-1, -1, 0, 0, 0, 0]
    inductance_xmax = [1, 1, 2π, 1, 1, 2π]
    
    val1, err = hcubature(
        inductance_cubature_func_1, 
        inductance_xmin,
        inductance_xmax;
        atol=abstol,
        rtol=reltol
    )
    val2, err = hcubature(
        inductance_cubature_func_2, 
        inductance_xmin,
        inductance_xmax;
        atol=abstol,
        rtol=reltol
    )
    val3, err = hcubature(
        inductance_cubature_func_3, 
        inductance_xmin,
        inductance_xmax;
        atol=abstol,
        rtol=reltol
    )
    val4, err = hcubature(
        inductance_cubature_func_4, 
        inductance_xmin,
        inductance_xmax;
        atol=abstol,
        rtol=reltol
    )
    val = val1 + val2 + val3 + val4
    sqrt_g_factor = 1 / 4
    prefactor = μ0 / (4π) * sqrt_g_factor * sqrt_g_factor
    return prefactor * val
end