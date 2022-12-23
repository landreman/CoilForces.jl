using Cubature

struct Coil
    curve::Curve
    current::Float64
    aminor::Float64
end

function analytic_force_per_unit_length(coil::Coil)
    # Assert that curve type is a CurveCircle:
    coil.curve::CurveCircle
    I = coil.current
    R = coil.curve.R0
    a = coil.aminor
    return μ0 * I * I / (4π * R) * (log(8 * R / a) - 0.75)
end

Biot_savart_prefactor = μ0 / (4π)

function d_B_d_ϕ(coil::Coil, ϕ, r_eval, regularization=0.0)
    Δr = r_eval - γ(coil.curve, ϕ)
    temp = normsq(Δr) + regularization
    denominator = temp * sqrt(temp)
    return coil.current * Biot_savart_prefactor * cross(dγdϕ(coil.curve, ϕ), Δr) / denominator
end

"""
    B_filament(coil::Coil, r_eval, nϕ, regularization=0.0)

Evaluate the Biot-Savart law for a coil in the approximation that the coil is an
infinitesmally thin filament.
"""
function B_filament_fixed(coil::Coil, r_eval, nϕ; regularization=0.0)
    dϕ = 2π / nϕ
    B = [0.0, 0.0, 0.0]
    for j in 1:nϕ
        ϕ = (j - 1) * dϕ
        B += d_B_d_ϕ(coil, ϕ, r_eval, regularization)
    end
    B *= dϕ
    return B
end

"""
Adaptive quadrature.
"""
function B_filament_adaptive(coil::Coil, r_eval; regularization=0.0, reltol=1e-8, abstol=1e-14)
    function Biot_Savart_integrand!(ϕ0, v) 
        v[:] = d_B_d_ϕ(coil, ϕ0, r_eval, regularization)
    end
    val, err = hquadrature(3, Biot_Savart_integrand!, 0, 2π, reltol=reltol, abstol=abstol)
    return val
end