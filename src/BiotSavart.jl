struct Coil
    curve::Curve
    current::Float64
    aminor::Float64
end

Biot_savart_prefactor = μ0 / (4π)

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
    B_filament_fixed(coil::Coil, r_eval, nϕ; regularization=0.0, drop_first_point=false)

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
    B_filament_adaptive(coil::Coil, r_eval; regularization=0.0, reltol=1e-8, abstol=1e-14)

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

function singularity_term(coil::Coil, ϕ)
    data = γ_and_2_derivatives(coil.curve, ϕ)
    r_prime = data[:, 2]
    dϕdℓ = 1 / norm(r_prime)
    δ = coil.aminor * coil.aminor / sqrt(ℯ)
    return (μ0 * coil.current / (4π) * dϕdℓ * dϕdℓ * dϕdℓ * (1 + log(sqrt(δ) * dϕdℓ / 8)) 
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
    δ = coil.aminor * coil.aminor / sqrt(ℯ)
    for j in 1:nϕ
        ϕ = (j - 1) * dϕ
        B += d_B_d_ϕ_singularity_subtracted(coil, ϕ, r_eval, δ, ϕ0, r_prime_prime_cross_r_prime, dℓdϕ_squared)
    end
    B = B * dϕ + singularity_term(coil, ϕ0)
    return B
end

"""
For a finite-thickness coil of arbitrary smooth shape, 
compute the integrand for evaluating B.

See 20221016-01 Numerical evaluation of B for finite thickness coil.lyx
"""
function B_finite_thickness_integrand(coil::Coil, ρ, θ, ϕ, r_eval)
    r = ρ * coil.aminor
    dℓdϕ, κ, τ, dr, tangent, normal, binormal = Frenet_frame(coil.curve, ϕ)
    cosθ = cos(θ)
    @. dr += (r * cosθ) * normal + (r * sin(θ)) * binormal - r_eval
    temp = 1 / (normsq(dr) + 1e-200)
    sqrtg = (1 - κ * r * cosθ) * ρ * dℓdϕ
    return (sqrtg * temp * sqrt(temp)) * cross(dr, tangent)
end

"""
Compute the magnetic field vector at a point with specified Cartesian
coordinates. In this version of the function, the prefactor μ0 I / (4 π^2) is
not included!
"""
function B_finite_thickness_normalized(coil::Coil, r_eval; reltol=1e-3, abstol=1e-5)
    function Biot_savart_cubature_func(xp)
        return B_finite_thickness_integrand(coil, xp[1], xp[2], xp[3], r_eval)
    end

    Biot_savart_xmin = [0, 0, 0]
    Biot_savart_xmax = [1, 2π, 2π]
    #Biot_savart_xmin = [0, -π, -π]
    #Biot_savart_xmax = [a, π, π]
    #Biot_savart_xmin = [0, 0.1 - π, 0.1 - π]
    #Biot_savart_xmax = [a, 0.1 + π, 0.1 + π]

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
function B_finite_thickness(coil::Coil, r_eval; reltol=1e-3, abstol=1e-5)
    prefactor = coil.current / (π) * Biot_savart_prefactor
    return prefactor * B_finite_thickness_normalized(coil, r_eval; reltol=reltol, abstol=abstol)
end