"""
Compute the force-per-unit-length for a finite-thickness coil.

ϕ: Curve parameter at which the force-per-unit-length will be computed.
"""
function force_finite_thickness(coil::Coil, ϕ; reltol=1e-3, abstol=1e-5)
    r0 = γ(coil.curve, ϕ)
    dℓdϕ, κ, τ, γ0, tangent, normal, binormal = Frenet_frame(coil.curve, ϕ)

    prefactor = coil.current / (π * coil.aminor * coil.aminor)

    function force_cubature_func!(xp, v)
        ρ = xp[1]
        θ = xp[2]
        cosθ = cos(θ)
        sqrtg = (1 - κ * ρ * cosθ) * ρ
        r_eval = r0 + ρ * cosθ * normal + ρ * sin(θ) * binormal
        B = B_finite_thickness(coil, r_eval, reltol=reltol, abstol=abstol)
        v[:] = sqrtg * cross(tangent, B)
    end

    force_xmin = [0, 0]
    force_xmax = [coil.aminor, 2π]
    
    val, err = hcubature(
        3,
        force_cubature_func!, 
        force_xmin,
        force_xmax,
        abstol=abstol,
        reltol=reltol)
    return prefactor * val
end