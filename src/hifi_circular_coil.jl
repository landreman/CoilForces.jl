"""
Tailored functions specifically for high-fidelity calculations with a circular
coil.
    
Ideas:
* Only compute the z component of B, since other components don't end up being
  needed.
* Do not include factors of μ0 and I in the integrands, to try to keep things O(1).
"""

"""
Returns just the z component of the Biot-Savart integrand

(x, 0, z) is the evaluation point.
(r, θ, ϕ) is the source location.
"""
function hifi_circular_coil_Biot_savart_z_integrand(R0, a, x, z, ρ, θ, ϕ)
    r = ρ * a
    cosϕ = cos(ϕ)
    sinϕ = sin(ϕ)
    R = R0 + r * cos(θ)
    sqrtg = ρ * R
    xp = R * cosϕ
    yp = R * sinϕ
    zp = r * sin(θ)
    dx = x - xp
    #dy = y - yp
    dy = -yp  # Since y=0
    dz = z - zp
    # Add a tiny number before dividing, so avoid divide-by-0 errors.
    dr_inv = 1 / (sqrt(dx * dx + dy * dy + dz * dz) + 1.0e-100)
    factor = dr_inv * dr_inv * dr_inv * sqrtg
    return -factor * (sinϕ * dy + cosϕ * dx)
end

"""
Compute B_z at a point with specified Cartesian coordinates. The prefactor 
I * μ0 / (4 π^2) is not included!
"""
function hifi_circular_coil_compute_Bz(R0, a, x, z; reltol=1e-3, abstol=1e-5)
    function Biot_savart_cubature_func(xp)
        return hifi_circular_coil_Biot_savart_z_integrand(R0, a, x, z, xp[1], xp[2], xp[3])
    end

    val, err = hcubature(
        Biot_savart_cubature_func, 
        [0, 0, 0],  # Lower integration bounds for (ρ, θ, ϕ)
        [1, 2π, 2π];  # Upper integration bounds for (ρ, θ, ϕ)
        atol=abstol,
        rtol=reltol,
        )
    return val
end

"""
Compute the x component of the force per unit length on a circular coil, using
the high fidelity model.
"""
function hifi_circular_coil_force(R0, a, I; reltol=1e-3, abstol=1e-5)
    """
    Given (r, θ) in the ϕ=0 plane, evaluate the integrand for computing the x
    component of the force, F_x:
    """
    function force_integrand(xx)
        ρ = xx[1]
        θ = xx[2]
        r = ρ * a
        R = R0 + r * cos(θ)
        # Evaluate B_z at (x, y, z)=(R, 0, r * sin(θ)):
        Bz = hifi_circular_coil_compute_Bz(R0, a, R, r * sin(θ); abstol=abstol, reltol=reltol)
        return ρ * (R / R0) * Bz
    end

    @time force_without_prefactors, force_err = hcubature(
        force_integrand, 
        [0, 0],  # Lower integration bounds for (ρ, θ)
        [1, 2π];  # Upper integration bounds for (ρ, θ)
        atol=abstol,
        rtol=reltol,
        )

    Biot_savart_prefactor = I * μ0 / (4 * π^2)
    force_prefactor = I / π
    return force_prefactor * Biot_savart_prefactor * force_without_prefactors
end