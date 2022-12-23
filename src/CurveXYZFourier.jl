struct CurveXYZFourier <: Curve
    xc::Vector{Float64}
    xs::Vector{Float64}
    yc::Vector{Float64}
    ys::Vector{Float64}
    zc::Vector{Float64}
    zs::Vector{Float64}
    n::Int
end

function CurveXYZFourier(xc, xs, yc, ys, zc, zs) 
    @assert length(xs) == length(xc)
    @assert length(yc) == length(xc)
    @assert length(ys) == length(xc)
    @assert length(zc) == length(xc)
    @assert length(zs) == length(xc)
    return CurveXYZFourier(xc, xs, yc, ys, zc, zs, length(xc))
end

function γ(c::CurveXYZFourier, ϕ)
    γ = [0.0, 0.0, 0.0]
    for m in 0 : (c.n - 1)
        cosfac = cos(m * ϕ)
        sinfac = sin(m * ϕ)

        γ += [
            c.xc[j] * cosfac + c.xs[j] * sinfac,
            c.yc[j] * cosfac + c.ys[j] * sinfac,
            c.zc[j] * cosfac + c.zs[j] * sinfac
        ]
    end
    return γ
end

function dγdϕ(c::CurveXYZFourier, ϕ)
    dγdϕ = [0.0, 0.0, 0.0]
    for m in 0 : (c.n - 1)
        cosfac = -m * sin(m * ϕ)
        sinfac = m * cos(m * ϕ)

        dγdϕ += [
            c.xc[j] * cosfac + c.xs[j] * sinfac,
            c.yc[j] * cosfac + c.ys[j] * sinfac,
            c.zc[j] * cosfac + c.zs[j] * sinfac
        ]
    end
    return dγdϕ
end

function d2γdϕ2(c::CurveXYZFourier, ϕ)
    d2γdϕ2 = [0.0, 0.0, 0.0]
    for m in 0 : (c.n - 1)
        cosfac = -m * m * cos(m * ϕ)
        sinfac = -m * m * sin(m * ϕ)

        d2γdϕ2 += [
            c.xc[j] * cosfac + c.xs[j] * sinfac,
            c.yc[j] * cosfac + c.ys[j] * sinfac,
            c.zc[j] * cosfac + c.zs[j] * sinfac
        ]
    end
    return d2γdϕ2
end


function d3γdϕ(c::CurveXYZFourier, ϕ)
    d3γdϕ3 = [0.0, 0.0, 0.0]
    for m in 0 : (c.n - 1)
        cosfac = m * m * m * sin(m * ϕ)
        sinfac = -m * m * m * cos(m * ϕ)

        d3γdϕ3 += [
            c.xc[j] * cosfac + c.xs[j] * sinfac,
            c.yc[j] * cosfac + c.ys[j] * sinfac,
            c.zc[j] * cosfac + c.zs[j] * sinfac
        ]
    end
    return d3γdϕ3
end
3