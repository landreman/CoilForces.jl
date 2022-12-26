abstract type Curve end

function tangent(c::Curve, ϕ)
    r_prime = dγdϕ(c, ϕ)
    return r_prime / norm(r_prime)    
end

function Frenet_frame(c::Curve, ϕ)
    # See
    # https://en.wikipedia.org/wiki/Frenet%E2%80%93Serret_formulas#Other_expressions_of_the_frame
    
    r_prime = dγdϕ(c, ϕ)
    r_prime_prime = d2γdϕ2(c, ϕ)
    r_prime_prime_prime = d3γdϕ3(c, ϕ)

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

    return differential_arclength, curvature, torsion, tangent, normal, binormal
end


function curvature(c::Curve, ϕ)
    r_prime = dγdϕ(c, γ)
    r_prime_prime = d2γdϕ2(c, γ)
    r_prime_prime_prime = d3γdϕ3(c, γ)

    norm_r_prime = norm(r_prime)
    differential_arclength = norm_r_prime
    r_prime_cross_r_prime_prime = cross(r_prime, r_prime_prime)
    norm_r_prime_cross_r_prime_prime = norm(r_prime_cross_r_prime_prime)
    curvature = (norm_r_prime_cross_r_prime_prime 
                / (norm_r_prime * norm_r_prime * norm_r_prime))
    return differential_arclength, curvature
end