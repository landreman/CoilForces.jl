abstract type Curve end

"""
Returns a tuple (position vector, tangent vector)
"""
function position_and_tangent(c::Curve, ϕ)
    data = γ_and_derivative(c, ϕ)
    r_prime = data[:, 2]
    return data[:, 1], r_prime / norm(r_prime)
end

function Frenet_frame(c::Curve, ϕ)
    # See
    # https://en.wikipedia.org/wiki/Frenet%E2%80%93Serret_formulas#Other_expressions_of_the_frame
    
    data = γ_and_3_derivatives(c, ϕ)
    r_prime = @view data[:, 2]
    r_prime_prime = @view data[:, 3]
    r_prime_prime_prime = @view data[:, 4]

    norm_r_prime = norm(r_prime)
    differential_arclength = norm_r_prime

    tangent = r_prime / norm_r_prime

    r_prime_cross_r_prime_prime = cross(r_prime, r_prime_prime)
    norm_r_prime_cross_r_prime_prime = norm(r_prime_cross_r_prime_prime)

    curvature = (norm_r_prime_cross_r_prime_prime 
                / (norm_r_prime * norm_r_prime * norm_r_prime))

    normal = cross(r_prime_cross_r_prime_prime, r_prime) / (norm_r_prime * norm_r_prime_cross_r_prime_prime)

    torsion = (dot(r_prime_cross_r_prime_prime, r_prime_prime_prime) 
                / (norm_r_prime_cross_r_prime_prime * norm_r_prime_cross_r_prime_prime))

    binormal = r_prime_cross_r_prime_prime / norm_r_prime_cross_r_prime_prime

    return differential_arclength, curvature, torsion, data[:, 1], tangent, normal, binormal
end

function Frenet_frame_without_torsion(c::Curve, ϕ)
    # See
    # https://en.wikipedia.org/wiki/Frenet%E2%80%93Serret_formulas#Other_expressions_of_the_frame
    
    data = γ_and_3_derivatives(c, ϕ)
    r_prime = @view data[:, 2]
    r_prime_prime = @view data[:, 3]
    r_prime_prime_prime = @view data[:, 4]

    norm_r_prime = norm(r_prime)
    differential_arclength = norm_r_prime

    tangent = r_prime / norm_r_prime

    r_prime_cross_r_prime_prime = cross(r_prime, r_prime_prime)
    norm_r_prime_cross_r_prime_prime = norm(r_prime_cross_r_prime_prime)

    curvature = (norm_r_prime_cross_r_prime_prime 
                / (norm_r_prime * norm_r_prime * norm_r_prime))

    normal = cross(r_prime_cross_r_prime_prime, r_prime) / (norm_r_prime * norm_r_prime_cross_r_prime_prime)

    return differential_arclength, curvature, data[:, 1], tangent, normal
end


function curvature(c::Curve, ϕ)
    data = γ_and_2_derivatives(c, ϕ)
    r_prime = @view data[:, 2]
    r_prime_prime = @view data[:, 3]

    norm_r_prime = norm(r_prime)
    differential_arclength = norm_r_prime
    r_prime_cross_r_prime_prime = cross(r_prime, r_prime_prime)
    norm_r_prime_cross_r_prime_prime = norm(r_prime_cross_r_prime_prime)
    curvature = (norm_r_prime_cross_r_prime_prime 
                / (norm_r_prime * norm_r_prime * norm_r_prime))
    return differential_arclength, curvature
end

function curve_length(c::Curve)
    function length_integrand(ϕ)
        data = γ_and_derivative(c, ϕ)
        r_prime = @view data[:, 2]
        return norm(r_prime)
    end

    val, err = hquadrature(
        length_integrand, 
        0,
        2π;
        atol=1e-12,
        rtol=1e-12
    )
    return val
end

function centroid(c::Curve)
    n = 1000
    x = zeros(3)
    denominator = 0.0
    for j in 1:n
        ϕ = (j - 1) * 2π / n
        data = γ_and_derivative(c, ϕ)
        d_l_d_ϕ = norm(data[:, 2])
        denominator += d_l_d_ϕ
        x += d_l_d_ϕ * data[:, 1]
    end
    return (1 / denominator) * x
end