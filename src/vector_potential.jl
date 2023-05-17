"""
For a finite-thickness coil of arbitrary smooth shape, 
compute the integrand for evaluating B.

See 20221016-01 Numerical evaluation of B for finite thickness coil.lyx
"""
function A_finite_thickness_integrand(coil::Coil, ρ, θ, ϕ, r_eval, regularization=1e-100)
    r = ρ * coil.aminor
    dℓdϕ, κ, τ, dr, tangent, normal, binormal = Frenet_frame(coil.curve, ϕ)
    sinθ, cosθ = sincos(θ)
    @. dr += (r * cosθ) * normal + (r * sinθ) * binormal - r_eval
    sqrtg = (1 - κ * r * cosθ) * ρ * dℓdϕ
    return (sqrtg * sqrt(1 / (normsq(dr) + regularization))) * tangent
end

"""
Compute the vector potential at a point with specified Cartesian
coordinates. In this version of the function, the prefactor μ0 I / (4 π^2) is
not included!
"""
function A_finite_thickness_normalized(coil::Coil, r_eval; reltol=1e-3, abstol=1e-5, ϕ_shift=0.0, θ_shift=0.0)
    function A_cubature_func(xp)
        return A_finite_thickness_integrand(coil, xp[1], xp[2], xp[3], r_eval)
    end

    Biot_savart_xmin = [0, θ_shift, ϕ_shift]
    Biot_savart_xmax = [1, θ_shift + 2π, ϕ_shift + 2π]

    val, err = hcubature(
        A_cubature_func, 
        Biot_savart_xmin,
        Biot_savart_xmax;
        atol=abstol,
        rtol=reltol,
    )
    return val
end

"""
Compute the vector potential at a point with specified Cartesian
coordinates. In this version of the function, the prefactor μ0 I / (4 π^2) is included.
"""
function A_finite_thickness(coil::Coil, r_eval; reltol=1e-3, abstol=1e-5, ϕ_shift=0.0, θ_shift=0.0)
    prefactor = μ0 * coil.current / (4 * π * π)
    return prefactor * A_finite_thickness_normalized(
        coil,
        r_eval;
        reltol,
        abstol,
        ϕ_shift,
        θ_shift,
    )
end
