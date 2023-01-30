"""
A curve which is just a circle in the z=0 plane, with radius R0.
"""
struct CurveCircle <: Curve
    R0::Float64
    buffer4::Matrix{Float64}
end

"""
Constructor.
"""
function CurveCircle(R0) 
    return CurveCircle(R0, zeros(3, 4))
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
    """
    return [
        (c.R0 * cosϕ) (-c.R0 * sinϕ) (-c.R0 * cosϕ) (c.R0 * sinϕ);
        (c.R0 * sinϕ) (c.R0 * cosϕ) (-c.R0 * sinϕ) (-c.R0 * cosϕ);
        0. 0. 0. 0.
    ]
    """
    """
    c.buffer4[:, :] .= [
        (c.R0 * cosϕ) (-c.R0 * sinϕ) (-c.R0 * cosϕ) (c.R0 * sinϕ);
        (c.R0 * sinϕ) (c.R0 * cosϕ) (-c.R0 * sinϕ) (-c.R0 * cosϕ);
        0. 0. 0. 0.
    ]
    """
    c.buffer4[1, 1] = c.R0 * cosϕ
    c.buffer4[1, 2] = -c.R0 * sinϕ
    c.buffer4[1, 3] = -c.R0 * cosϕ
    c.buffer4[1, 4] = c.R0 * sinϕ
    c.buffer4[2, 1] = c.R0 * sinϕ
    c.buffer4[2, 2] = c.R0 * cosϕ
    c.buffer4[2, 3] = -c.R0 * sinϕ
    c.buffer4[2, 4] = -c.R0 * cosϕ
    c.buffer4[3, :] .= 0.0
    return c.buffer4
end
