struct CurveXYZFourier <: Curve
    cos_coeffs::Matrix{Float64}
    sin_coeffs::Matrix{Float64}
    n::Int
    buffer_33::Matrix{Float64}
    buffer_34::Matrix{Float64}
    cos_m_factors::Matrix{Float64}
    sin_m_factors::Matrix{Float64}
end

"""
Constructor, which checks the lengths of the input vectors of Fourier amplitudes.
"""
function CurveXYZFourier(xc, xs, yc, ys, zc, zs) 
    @assert length(xs) == length(xc)
    @assert length(yc) == length(xc)
    @assert length(ys) == length(xc)
    @assert length(zc) == length(xc)
    @assert length(zs) == length(xc)
    cos_coeffs = [
        xc';
        yc';
        zc';
    ]
    sin_coeffs = [
        xs';
        ys';
        zs';
    ]
    n = length(xc)
    return CurveXYZFourier(
        cos_coeffs, 
        sin_coeffs, 
        n, 
        zeros(3, 3),
        zeros(3, 4),
        zeros(n, 4),
        zeros(n, 4),
    )
end

function γ(c::CurveXYZFourier, ϕ)
    cos_m_factors = [cos(m * ϕ) for m in 0:(c.n - 1)]
    sin_m_factors = [sin(m * ϕ) for m in 0:(c.n - 1)]
    return c.cos_coeffs * cos_m_factors + c.sin_coeffs * sin_m_factors
end

"""
    γ_and_derivative(c::CurveXYZFourier, ϕ)

Returns a 3 × 4 matrix containing
[γ, dγ/dϕ, d^2γ/dϕ^2, d^3γ/dϕ^3]
"""
function γ_and_derivative(c::CurveXYZFourier, ϕ)

    cos_m_factors = ones(c.n, 2)
    sin_m_factors = ones(c.n, 2)

    for j in 1:c.n
        m = j - 1
        sin_mϕ, cos_mϕ = sincos(m * ϕ)

        cos_m_factors[j, 1] = cos_mϕ
        cos_m_factors[j, 2] = -m * sin_mϕ

        sin_m_factors[j, 1] = sin_mϕ
        sin_m_factors[j, 2] = m * cos_mϕ
    end

    data = c.cos_coeffs * cos_m_factors + c.sin_coeffs * sin_m_factors
    
    return data
end

"""
    γ_and_2_derivatives(c::CurveXYZFourier, ϕ)

Returns a 3 × 4 matrix containing
[γ, dγ/dϕ, d^2γ/dϕ^2, d^3γ/dϕ^3]
"""
function γ_and_2_derivatives(c::CurveXYZFourier, ϕ)

    cos_m_factors = ones(c.n, 3)
    sin_m_factors = ones(c.n, 3)

    for j in 1:c.n
        m = j - 1
        sin_mϕ, cos_mϕ = sincos(m * ϕ)

        cos_m_factors[j, 1] = cos_mϕ
        cos_m_factors[j, 2] = -m * sin_mϕ
        cos_m_factors[j, 3] = -m * m * cos_mϕ

        sin_m_factors[j, 1] = sin_mϕ
        sin_m_factors[j, 2] = m * cos_mϕ
        sin_m_factors[j, 3] = -m * m * sin_mϕ
    end

    data = c.cos_coeffs * cos_m_factors + c.sin_coeffs * sin_m_factors
    
    return data
end

"""
    γ_and_3_derivatives(c::CurveXYZFourier, ϕ)

Returns a 3 × 4 matrix containing
[γ, dγ/dϕ, d^2γ/dϕ^2, d^3γ/dϕ^3]
"""
function γ_and_3_derivatives(c::CurveXYZFourier, ϕ)

    for j in 1:c.n
        m = j - 1
        sin_mϕ, cos_mϕ = sincos(m * ϕ)

        c.cos_m_factors[j, 1] = cos_mϕ
        c.cos_m_factors[j, 2] = -m * sin_mϕ
        c.cos_m_factors[j, 3] = -m * m * cos_mϕ
        c.cos_m_factors[j, 4] = m * m * m * sin_mϕ

        c.sin_m_factors[j, 1] = sin_mϕ
        c.sin_m_factors[j, 2] = m * cos_mϕ
        c.sin_m_factors[j, 3] = -m * m * sin_mϕ
        c.sin_m_factors[j, 4] = -m * m * m * cos_mϕ
    end

    #c.buffer4 = c.cos_coeffs * cos_m_factors + c.sin_coeffs * sin_m_factors
    LinearAlgebra.mul!(c.buffer_34, c.cos_coeffs, c.cos_m_factors)
    LinearAlgebra.mul!(c.buffer_34, c.sin_coeffs, c.sin_m_factors, 1, 1)
    
    return c.buffer_34
end

"""
Given any Curve, returns a CurveXYZFourier that describes a circle (translated
and rotated) which is a best-fit approximation to the original curve at a
desired point ϕ.
"""
function fit_circle(curve::Curve, ϕ)
    _, curvature, _, position, tangent, normal, _ = Frenet_frame(curve, ϕ)
    sinϕ, cosϕ = sincos(ϕ)

    # For the formulas that follow, see
    # 20230303-01 High fidelity force calculation subtracting singularity of best fit circle.lyx
    return CurveXYZFourier(
        # xc:
        [position[1] + normal[1] / curvature, -(sinϕ * tangent[1] + cosϕ * normal[1]) / curvature],
        # xs:
        [0, (cosϕ * tangent[1] - sinϕ * normal[1]) / curvature],
        # yc:
        [position[2] + normal[2] / curvature, -(sinϕ * tangent[2] + cosϕ * normal[2]) / curvature],
        # ys:
        [0, (cosϕ * tangent[2] - sinϕ * normal[2]) / curvature],
        # zc:
        [position[3] + normal[3] / curvature, -(sinϕ * tangent[3] + cosϕ * normal[3]) / curvature],
        # zs:
        [0, (cosϕ * tangent[3] - sinϕ * normal[3]) / curvature],
    )
end