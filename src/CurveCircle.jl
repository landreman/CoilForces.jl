"""
A curve which is just a circle in the z=0 plane, with radius R0.
"""
struct CurveCircle <: Curve
    R0::Float64
end

function γ(c::CurveCircle, ϕ)
    return [c.R0 * cos(ϕ), c.R0 * sin(ϕ), 0.0]
end

function dγdϕ(c::CurveCircle, ϕ)
    return [-c.R0 * sin(ϕ), c.R0 * cos(ϕ), 0.0]
end

function d2γdϕ2(c::CurveCircle, ϕ)
    return [-c.R0 * cos(ϕ), -c.R0 * sin(ϕ), 0.0]
end

function d3γdϕ3(c::CurveCircle, ϕ)
    return [c.R0 * sin(ϕ), -c.R0 * cos(ϕ), 0.0]
end
