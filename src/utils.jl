const μ0 = 4π * (1e-7)

dot(v1, v2) = v1[1] * v2[1] + v1[2] * v2[2] + v1[3] * v2[3]

normsq(v) = v[1] * v[1] + v[2] * v[2] + v[3] * v[3]

norm(v) = sqrt(normsq(v))

function cross(v1, v2)
    return [v1[2] * v2[3] - v1[3] * v2[2],
            v1[3] * v2[1] - v1[1] * v2[3],
            v1[1] * v2[2] - v1[2] * v2[1]]
end
