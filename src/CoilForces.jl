module CoilForces

export μ0, dot, cross
export CurveCircle, CurveXYZFourier
export γ, dγdϕ, d2γdϕ2, d3γdϕ3, curvature_torsion
export Coil
export analytic_force_per_unit_length

include("utils.jl")
include("Curve.jl")
include("CurveCircle.jl")
include("CurveXYZFourier.jl")
include("BiotSavart.jl")

end
