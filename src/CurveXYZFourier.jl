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
        cos_mϕ = cos(m * ϕ)
        sin_mϕ = sin(m * ϕ)

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
        cos_mϕ = cos(m * ϕ)
        sin_mϕ = sin(m * ϕ)

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
        cos_mϕ = cos(m * ϕ)
        sin_mϕ = sin(m * ϕ)

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
