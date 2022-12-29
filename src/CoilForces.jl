module CoilForces

using Cubature
using Plots
using Printf

export μ0, dot, cross, norm, normsq
export CurveCircle, CurveXYZFourier
export γ, dγdϕ, d2γdϕ2, d3γdϕ3, tangent, Frenet_frame
export Coil
export d_B_d_ϕ, B_filament_fixed, B_filament_adaptive, analytic_force_per_unit_length
export d_B_d_ϕ_singularity_subtracted, B_singularity_subtraction_fixed
export B_finite_thickness
export get_curve
export plot_force_for_HSX, plot_integrand
export force_finite_thickness

include("utils.jl")
include("Curve.jl")
include("CurveCircle.jl")
include("CurveXYZFourier.jl")
include("configs.jl")
include("BiotSavart.jl")
include("results.jl")
include("force.jl")

end
