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
There are two types for representing the winding pack angle: FrameCircle and FrameCentroid.
These are both subtypes of an abstract type Frame.
"""
abstract type Frame end

struct FrameCircle <: Frame
end

function get_frame(frame::FrameCircle, ϕ)
    sinϕ, cosϕ = sincos(ϕ)
    return [-cosϕ, -sinϕ, 0.0], [0.0, 0.0, 1.0]
end

struct FrameCentroid <: Frame
    centroid::Vector{Float64}
end

"""
Constructor, which takes a curve
"""
function FrameCentroid(curve::Curve)
    return FrameCentroid(centroid(curve))
end

function get_frame(frame::FrameCentroid, ϕ)
    # Replace this content with the real function eventually
    sinϕ, cosϕ = sincos(ϕ)
    return [-cosϕ, -sinϕ, 0.0], [0.0, 0.0, 1.0]
end

struct FrameRotated <: Frame
    frame::Frame
    cos_angle::Float64
    sin_angle::Float64
end

function FrameRotated(frame::Frame, angle)
    return FrameRotated(frame, cos(angle), sin(angle))
end

function get_frame(frame::FrameRotated, ϕ)
    p, q = get_frame(frame.frame, ϕ)
    return (
        p * frame.cos_angle + q * frame.sin_angle,
        -p * frame.sin_angle, q * frame_cos_angle
    )
end

struct CoilRectangularXSection <: Coil
    curve::Curve
    current::Float64
    a::Float64
    b::Float64
    frame::Frame
end

function get_κ1_κ2(p, q, normal, κ)
    cos_α = dot(p, normal)
    minus_sin_α = dot(q, normal)
    κ1 = κ * cos_α
    κ2 = κ * minus_sin_α
    return κ1, κ2
end