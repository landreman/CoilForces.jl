Biot_savart_prefactor = μ0 / (4π)

sqrt_1_over_e = 1 / √ℯ

function compute_regularization(coil::CoilCircularXSection)
    return coil.aminor * coil.aminor * sqrt_1_over_e
end

function compute_regularization(coil::CoilRectangularXSection)
    return coil.a * coil.b * CoilForces.rectangular_xsection_δ(coil.a, coil.b)
end

function d_B_d_ϕ(coil::Coil, ϕ, r_eval; regularization=0.0)
    data = γ_and_derivative(coil.curve, ϕ)
    Δr = r_eval - data[:, 1]
    temp = normsq(Δr) + regularization
    denominator = temp * sqrt(temp)
    return coil.current * Biot_savart_prefactor * cross(data[:, 2], Δr) / denominator
end

"""
ϕ: curve parameter for the incremental current that contributes to B.
ϕ0: curve parameter at which we are evaluating B.
"""
function d_B_d_ϕ_singularity_subtracted(coil::Coil, ϕ, r_eval, regularization, ϕ0, r_prime_prime_cross_r_prime, dℓdϕ_squared)
    data = γ_and_derivative(coil.curve, ϕ)
    Δr = r_eval - data[:, 1]
    temp = normsq(Δr) + regularization
    denominator = temp * sqrt(temp)
    cos_fac = 2 - 2 * cos(ϕ - ϕ0)
    temp2 = cos_fac * dℓdϕ_squared + regularization
    denominator2 = temp2 * sqrt(temp2)
    return coil.current * Biot_savart_prefactor * (
        cross(data[:, 2], Δr) / denominator
        + (0.5 * cos_fac / denominator2) * r_prime_prime_cross_r_prime
    )
end


"""
    B_filament_fixed(coil::CoilCircularXSection, r_eval, nϕ; regularization=0.0, drop_first_point=false)

Evaluate the Biot-Savart law for a coil in the approximation that the coil is an
infinitesmally thin filament. Use quadrature on a fixed uniform grid with
specified number of points, nϕ.
"""
function B_filament_fixed(coil::Coil, r_eval, nϕ; regularization=0.0, drop_first_point=false)
    dϕ = 2π / nϕ
    B = [0.0, 0.0, 0.0]
    if drop_first_point
        first_point = 2
    else
        first_point = 1
    end
    for j in first_point:nϕ
        ϕ = (j - 1) * dϕ
        B += d_B_d_ϕ(coil, ϕ, r_eval, regularization=regularization)
    end
    B *= dϕ
    return B
end

"""
    B_filament_adaptive(coil::CoilCircularXSection, r_eval; regularization=0.0, reltol=1e-8, abstol=1e-14)

Evaluate the Biot-Savart law for a coil in the approximation that the coil is an
infinitesmally thin filament. Use adaptive quadrature.
"""
function B_filament_adaptive(coil::Coil, r_eval; regularization=0.0, reltol=1e-8, abstol=1e-14)
    function Biot_Savart_integrand(ϕ0) 
        return d_B_d_ϕ(coil, ϕ0, r_eval, regularization=regularization)
    end
    val, err = hquadrature(Biot_Savart_integrand, 0, 2π; rtol=reltol, atol=abstol)
    return val
end

"""
This function computes the term to be added back in the singularity-subtraction
trick for regularized Biot-Savart, i.e. the piece that was integrated analytically.
"""
function singularity_term(coil::Coil, ϕ)
    data = γ_and_2_derivatives(coil.curve, ϕ)
    r_prime = data[:, 2]
    dϕdℓ = 1 / norm(r_prime)
    return (μ0 * coil.current / (8π) * dϕdℓ * dϕdℓ * dϕdℓ 
        * (2 + log(compute_regularization(coil) * dϕdℓ * dϕdℓ / 64)) 
        * cross(data[:, 3], r_prime))
end

"""
    B_singularity_subtraction_fixed(coil::Coil, ϕ0, nϕ)

Evaluate the Biot-Savart law for a coil in the approximation that the coil is an
infinitesmally thin filament. Use quadrature on a fixed uniform grid with
specified number of points, nϕ.

ϕ0: curve parameter at which to evaluate B.
"""
function B_singularity_subtraction_fixed(coil::Coil, ϕ0, nϕ)
    dϕ = 2π / nϕ
    B = [0.0, 0.0, 0.0]
    
    data = γ_and_2_derivatives(coil.curve, ϕ0)
    r_eval = data[:, 1]
    r_prime = data[:, 2]
    r_prime_prime = data[:, 3]
    
    dℓdϕ_squared = normsq(r_prime)
    r_prime_prime_cross_r_prime = cross(r_prime_prime, r_prime)
    regularization_to_use = compute_regularization(coil)
    for j in 1:nϕ
        ϕ = (j - 1) * dϕ
        B += d_B_d_ϕ_singularity_subtracted(coil, ϕ, r_eval, regularization_to_use, ϕ0, r_prime_prime_cross_r_prime, dℓdϕ_squared)
    end
    B = B * dϕ + singularity_term(coil, ϕ0)
    return B
end

"""
For a finite-thickness coil of arbitrary smooth shape, 
compute the integrand for evaluating B.

See 20221016-01 Numerical evaluation of B for finite thickness coil.lyx
"""
function B_finite_thickness_integrand(coil::CoilCircularXSection, ρ, θ, ϕ, r_eval, regularization=1e-100)
    r = ρ * coil.aminor
    dℓdϕ, κ, τ, dr, tangent, normal, binormal = Frenet_frame(coil.curve, ϕ)
    sinθ, cosθ = sincos(θ)
    @. dr += (r * cosθ) * normal + (r * sinθ) * binormal - r_eval
    #temp = 1 / (normsq(dr) + 1e-10)
    temp = 1 / (normsq(dr) + regularization)
    #temp = 1 / (normsq(dr) + 1e-200)
    sqrtg = (1 - κ * r * cosθ) * ρ * dℓdϕ
    return (sqrtg * temp * sqrt(temp)) * cross(dr, tangent)
end

"""
For a finite-thickness coil of arbitrary smooth shape, 
compute the integrand for evaluating B.
"""
function B_finite_thickness_integrand(coil::CoilRectangularXSection, u, v, ϕ, r_eval, regularization=1e-100)
    u_a_over_2 = 0.5 * u * coil.a
    v_b_over_2 = 0.5 * v * coil.b
    dℓdϕ, κ, dr, tangent, normal = Frenet_frame_without_torsion(coil.curve, ϕ)
    p, q = get_frame(coil.frame, ϕ, dr, tangent, normal)
    κ1, κ2 = CoilForces.get_κ1_κ2(p, q, normal, κ)
    # Note that sqrtg in the next line does not contain ab/4. We consider that
    # factor to be in front of the overall integral.
    sqrtg = (dℓdϕ 
        * (1 - u_a_over_2 * κ1 - v_b_over_2 * κ2))
    @. dr += u_a_over_2 * p + v_b_over_2 * q - r_eval
    temp = 1 / (normsq(dr) + regularization)
    return (sqrtg * temp * sqrt(temp)) * cross(dr, tangent)
end

"""
This version of the function takes cosθ and sinθ instead of θ. This improves
efficiency for singularity-subtraction calculations so the cos and sin do not
need to be recalculated.
"""
function B_finite_thickness_integrand_sincos(coil::CoilCircularXSection, ρ, cosθ, sinθ, ϕ, r_eval)
    r = ρ * coil.aminor
    dℓdϕ, κ, τ, dr, tangent, normal, binormal = Frenet_frame(coil.curve, ϕ)
    @. dr += (r * cosθ) * normal + (r * sinθ) * binormal - r_eval
    #temp = 1 / (normsq(dr) + 1e-30)
    temp = 1 / (normsq(dr) + 1e-100)
    #temp = 1 / (normsq(dr) + 1e-200)
    sqrtg = (1 - κ * r * cosθ) * ρ * dℓdϕ
    return (sqrtg * temp * sqrt(temp)) * cross(dr, tangent)
end

"""
Compute the magnetic field vector at a point with specified Cartesian
coordinates. In this version of the function, the prefactor μ0 I / (4 π^2) is
not included!
"""
function B_finite_thickness_normalized(coil::CoilCircularXSection, r_eval; reltol=1e-3, abstol=1e-5, ϕ_shift=0.0, θ_shift=0.0)
    function Biot_savart_cubature_func(xp)
        return B_finite_thickness_integrand(coil, xp[1], xp[2], xp[3], r_eval)
    end

    Biot_savart_xmin = [0, θ_shift, ϕ_shift]
    Biot_savart_xmax = [1, θ_shift + 2π, ϕ_shift + 2π]

    val, err = hcubature(
        Biot_savart_cubature_func, 
        Biot_savart_xmin,
        Biot_savart_xmax;
        atol=abstol,
        rtol=reltol,
        #maxevals=5000000,
        #initdiv=10,
    )
    #print("Number of function evals: ", myindex)
    return val
end

"""
Compute the magnetic field vector at a point with specified Cartesian
coordinates. In this version of the function, the prefactor μ0 I / (16 π) is
not included!
"""
function B_finite_thickness_normalized(coil::CoilRectangularXSection, r_eval; reltol=1e-3, abstol=1e-5, ϕ_shift=0.0)
    function Biot_savart_cubature_func(xp)
        return B_finite_thickness_integrand(coil, xp[1], xp[2], xp[3], r_eval)
    end

    Biot_savart_xmin = [-1, -1, ϕ_shift]
    Biot_savart_xmax = [1, 1, ϕ_shift + 2π]

    val, err = hcubature(
        Biot_savart_cubature_func, 
        Biot_savart_xmin,
        Biot_savart_xmax;
        atol=abstol,
        rtol=reltol,
    )
    return val
end

"""
Compute the magnetic field vector at a point with specified Cartesian
coordinates. In this version of the function, the prefactor μ0 I / (4 π^2) is
not included. Also, Siena's trick of subtracting the contribution from the
best-fit circular coil is used.
"""
function B_finite_thickness_singularity_subtraction(coil::CoilCircularXSection, best_fit_circular_coil::CoilCircularXSection, r_eval; reltol=1e-3, abstol=1e-5, ϕ_shift=0.0, θ_shift=0.0)
    function Biot_savart_cubature_func(xp)
        sinθ, cosθ = sincos(xp[2])
        return (
            B_finite_thickness_integrand_sincos(coil, xp[1], cosθ, sinθ, xp[3], r_eval)
            - B_finite_thickness_integrand_sincos(best_fit_circular_coil, xp[1], cosθ, sinθ, xp[3], r_eval)
        )
    end

    Biot_savart_xmin = [0, θ_shift, ϕ_shift]
    Biot_savart_xmax = [1, θ_shift + 2π, ϕ_shift + 2π]

    val, err = hcubature(
        Biot_savart_cubature_func, 
        Biot_savart_xmin,
        Biot_savart_xmax;
        atol=abstol,
        rtol=reltol)
    return val
end


"""
Compute the magnetic field vector at a point with specified Cartesian
coordinates. In this version of the function, the prefactor μ0 I / (4 π^2) is
included.
"""
function B_finite_thickness(coil::CoilCircularXSection, r_eval; reltol=1e-3, abstol=1e-5, ϕ_shift=0.0, θ_shift=0.0)
    prefactor = coil.current / (π) * Biot_savart_prefactor
    return prefactor * B_finite_thickness_normalized(coil, r_eval; reltol=reltol, abstol=abstol, ϕ_shift=ϕ_shift, θ_shift=θ_shift)
end

"""
Compute the magnetic field vector at a point with specified Cartesian
coordinates. In this version of the function, the prefactor μ0 I / (16 π) is
included.
"""
function B_finite_thickness(coil::CoilRectangularXSection, r_eval; reltol=1e-3, abstol=1e-5, ϕ_shift=0.0)
    prefactor = μ0 * coil.current / (16 * π)
    return prefactor * B_finite_thickness_normalized(coil, r_eval; reltol=reltol, abstol=abstol, ϕ_shift=ϕ_shift)
end

"""
Returns eq (124) in
20230326-01_B_in_conductor_for_a_noncircular_finite_thickness_coil.pdf
"""
function B_local(coil::CoilCircularXSection, curvature, normal, binormal, ρ, θ)
    return (
        μ0 * coil.current * ρ / (2π * coil.aminor) * (-normal * sin(θ) + binormal * cos(θ))
        + μ0 * coil.current * curvature / (8π) * (
            -0.5 * ρ^2 * sin(2θ) * normal
            + (1.5 + ρ^2 * (-1 + 0.5 * cos(2θ))) * binormal
        )
    )
end

function rectangular_xsection_G(p, q)
    return q * atan(p / q) + 0.5 * p * log(1 + q * q / (p * p))
end

"""
Returns the field of an infinite straight wire with rectangular cross-section
"""
function rectangular_xsection_B0(coil::CoilRectangularXSection, u, v)
    B0_q = 0
    B0_p = 0
    for su in [-1, 1]
        adu = coil.a * (u - su)
        for sv in [-1, 1]
            bdv = coil.b * (v - sv)
            B0_q += su * sv * rectangular_xsection_G(bdv, adu)
            B0_p -= su * sv * rectangular_xsection_G(adu, bdv)
        end
    end
    prefactor = μ0 * coil.current / (4π * coil.a * coil.b)
    return prefactor * B0_p, prefactor * B0_q
end

function rectangular_xsection_K(a, b, κ1, κ2, U, V)
    temp = a * U * U / b + b * V * V / a
    log_factor = log(temp)
    K_p = (
        -2 * U * V * (-κ2) * log_factor
        - κ1 * temp * log_factor
        + 4 * a * U * U * κ2 / b * atan(b * V / (a * U))
    )
    K_q = (
        -2 * U * V * κ1 * log_factor
        + κ2 * temp * log_factor
        - 4 * b * V * V * κ1 / a * atan(a * U / (b * V))
    )
    return K_p, K_q
end

function rectangular_xsection_Bκ(coil::CoilRectangularXSection, κ1, κ2, u, v)
    Bκ_p = 0
    Bκ_q = 0
    for su in [-1, 1]
        for sv in [-1, 1]
            K_p, K_q = rectangular_xsection_K(coil.a, coil.b, κ1, κ2, u - su, v - sv)
            Bκ_p += su * sv * K_p
            Bκ_q += su * sv * K_q
        end
    end
    prefactor = μ0 * coil.current / (64π)
    return prefactor * Bκ_p, prefactor * Bκ_q
end