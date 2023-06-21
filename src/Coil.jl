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

function get_frame(frame::FrameCircle, ϕ, position, tangent, normal)
    return normal, cross(tangent, normal)
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

function get_frame(frame::FrameCentroid, ϕ, position, tangent, normal)
    # Section 3.1 of Singh (2020)
    p = position - frame.centroid
    temp = dot(p, tangent)
    @. p -= temp * tangent
    temp2 = norm(p)
    p = (1 / temp2) * p
    return p, cross(tangent, p)
end

struct FrameRotated <: Frame
    frame::Frame
    cos_angle::Float64
    sin_angle::Float64
end

function FrameRotated(frame::Frame, angle)
    return FrameRotated(frame, cos(angle), sin(angle))
end

function get_frame(frame::FrameRotated, ϕ, position, tangent, normal)
    p, q = get_frame(frame.frame, ϕ, position, tangent, normal)
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

"""
For a coil with rectangular cross-section, save the shapes of the coil edges to
a CSV file. This is useful for plotting the coil shape with an outside package
such as mayavi.
"""
function save(coil::CoilRectangularXSection, filename, n)
    u = [1, 1, -1, -1]
    v = [1, -1, -1, 1]
    open(filename, "w") do file
        write(file, "x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4\n")
        for j in 1:n
            ϕ = 2π * j / n
            dℓdϕ, κ, rc, tangent, normal = Frenet_frame_without_torsion(coil.curve, ϕ)
            p, q = get_frame(coil.frame, ϕ, rc, tangent, normal)
            first_point = true
            for k in 1:4
                r = rc + u[k] * p * coil.a / 2 + v[k] * q * coil.b / 2
                if !first_point
                    write(file, ", ")
                end
                write(file, "$(r[1]), $(r[2]), $(r[3])")
                first_point = false
            end
            write(file, "\n")
        end
    end
end

function rectangular_xsection_k(a, b)
    return (
        -(a^4 - 6 * a^2 * b^2 + b^4) / (6 * a^2 * b^2) * log(a / b + b / a)
        + b * b / (6 * a * a) * log(b / a)
        + a * a / (6 * b * b) * log(a / b)
        + (4.0 * b) / (3 * a) * atan(a / b)
        + (4.0 * a) / (3 * b) * atan(b / a)
    )
end

function rectangular_xsection_δ(a, b)
    return exp(-(25.0 / 6) + rectangular_xsection_k(a, b))
end