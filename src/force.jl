"""
Compute the force-per-unit-length for a finite-thickness coil.

ϕ: Curve parameter at which the force-per-unit-length will be computed.
"""
function force_finite_thickness(coil::Coil, ϕ; reltol=1e-3, abstol=1e-5)
    r0 = γ(coil.curve, ϕ)
    dℓdϕ, κ, τ, γ0, tangent, normal, binormal = Frenet_frame(coil.curve, ϕ)

    function force_cubature_func(xp)
        ρ = xp[1]
        θ = xp[2]
        r = ρ * coil.aminor
        cosθ = cos(θ)
        sqrtg = (1 - κ * r * cosθ) * ρ
        r_eval = r0 + r * cosθ * normal + r * sin(θ) * binormal
        B = B_finite_thickness_normalized(coil, r_eval, reltol=reltol, abstol=abstol)
        return sqrtg * cross(tangent, B)
    end

    force_xmin = [0, 0]
    force_xmax = [1, 2π]
    
    val, err = hcubature(
        force_cubature_func, 
        force_xmin,
        force_xmax;
        atol=abstol,
        rtol=reltol
    )
    Biot_savart_prefactor = coil.current * μ0 / (4 * π^2)
    force_prefactor = coil.current / π
    return Biot_savart_prefactor * force_prefactor * val
end

"""
Compute the self-force per unit length, using the locally circular
approximation, Garren & Chen eq (34).
"""
function force_locally_circular_approximation(coil::Coil, ϕ)
    differential_arclength, curvature, torsion, position, tangent, normal, binormal = Frenet_frame(coil.curve, ϕ)
    circular_coil = Coil(CurveCircle(1 / curvature), coil.current, coil.aminor)
    return -normal * analytic_force_per_unit_length(circular_coil)
end