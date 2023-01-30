"""
A curve which is just a circle in the z=0 plane, with radius R0.
"""
struct CurveCircle <: Curve
    R0::Float64
end

function γ(c::CurveCircle, ϕ)
    return [c.R0 * cos(ϕ), c.R0 * sin(ϕ), 0.0]
end

"""
    γ_and_derivative(c::CurveXYZFourier, ϕ)

Returns a 3 × 2 matrix containing
[γ, dγ/dϕ]
"""
function γ_and_derivative(c::CurveCircle, ϕ)
    cosϕ = cos(ϕ)
    sinϕ = sin(ϕ)
    return [
        (c.R0 * cosϕ) (-c.R0 * sinϕ);
        (c.R0 * sinϕ) (c.R0 * cosϕ);
        0. 0.
    ]
end

"""
    γ_and_2_derivatives(c::CurveXYZFourier, ϕ)

Returns a 3 × 3 matrix containing
[γ, dγ/dϕ, d^2γ/dϕ^2]
"""
function γ_and_2_derivatives(c::CurveCircle, ϕ)
    cosϕ = cos(ϕ)
    sinϕ = sin(ϕ)
    return [
        (c.R0 * cosϕ) (-c.R0 * sinϕ) (-c.R0 * cosϕ);
        (c.R0 * sinϕ) (c.R0 * cosϕ) (-c.R0 * sinϕ);
        0. 0. 0.
    ]
end

"""
    γ_and_3_derivatives(c::CurveXYZFourier, ϕ)

Returns a 3 × 4 matrix containing
[γ, dγ/dϕ, d^2γ/dϕ^2, d^3γ/dϕ^3]
"""
function γ_and_3_derivatives(c::CurveCircle, ϕ)
    cosϕ = cos(ϕ)
    sinϕ = sin(ϕ)
    return [
        (c.R0 * cosϕ) (-c.R0 * sinϕ) (-c.R0 * cosϕ) (c.R0 * sinϕ);
        (c.R0 * sinϕ) (c.R0 * cosϕ) (-c.R0 * sinϕ) (-c.R0 * cosϕ);
        0. 0. 0. 0.
    ]
end
