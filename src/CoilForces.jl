module CoilForces

export μ0, dot, cross
export CurveCircle, CurveXYZFourier
export γ, dγdϕ, d2γdϕ2, d3γdϕ3, curvature_torsion
export Coil
export B_filament_fixed, B_filament_adaptive, analytic_force_per_unit_length
export get_curve

include("utils.jl")
include("Curve.jl")
include("CurveCircle.jl")
include("CurveXYZFourier.jl")
include("configs.jl")
include("BiotSavart.jl")

end
