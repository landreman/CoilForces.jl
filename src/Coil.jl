"""
There are two types of Coil: CoilCircularXSection and CoilRectangularXSection.
These are both subtypes of an abstract type Coil.
"""
abstract type Coil end

struct CoilCircularXSection <: Coil
    curve::Curve
    current::Float64
    aminor::Float64
end

"""
There are two types for representing the winding pack angle: WindingPackAngleZero and WindingPackAngleFourier.
These are both subtypes of an abstract type WindingPackAngle.
"""
abstract type WindingPackAngle end

struct WindingPackAngleZero <: WindingPackAngle
end

function get_winding_pack_angle(wpa::WindingPackAngleZero, ϕ)
    return 0.0
end

struct WindingPackAngleFourier <: WindingPackAngle
    cos_coeffs::Vector{Float64}
    sin_coeffs::Vector{Float64}
    n::Int
    cos_m_factors::Vector{Float64}
    sin_m_factors::Vector{Float64}
end

"""
Constructor, which checks the lengths of the input vectors of Fourier amplitudes.
"""
function WindingPackAngleFourier(cos_coeffs, sin_coeffs) 
    @assert length(cos_coeffs) == length(sin_coeffs)
    n = length(cos_coeffs)
    return WindingPackAngleFourier(
        cos_coeffs, 
        sin_coeffs, 
        n, 
        zeros(n),
        zeros(n),
    )
end

function get_winding_pack_angle(wpa::WindingPackAngleFourier, ϕ)
    for j in 1:wpa.n
        wpa.sin_m_factors[j], wpa.cos_m_factors[j] = sincos((j - 1) * ϕ)
    end
    return (
        LinearAlgebra.dot(wpa.cos_m_factors, wpa.cos_coeffs) 
        + LinearAlgebra.dot(wpa.sin_m_factors, wpa.sin_coeffs)
    )
end

struct CoilRectangularXSection <: Coil
    curve::Curve
    current::Float64
    a::Float64
    b::Float64
    winding_pack_angle::WindingPackAngle
end
