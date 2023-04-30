module CoilForces

import LinearAlgebra
using HCubature
using Plots
using Printf
using Dates
using Interpolations
using CSV

export μ0, dot, cross, norm, normsq
export CurveCircle, CurveXYZFourier, fit_circle
export γ, γ_and_derivative, γ_and_2_derivatives, γ_and_3_derivatives, tangent, Frenet_frame, curve_length
export Coil
export d_B_d_ϕ, B_filament_fixed, B_filament_adaptive
export d_B_d_ϕ_singularity_subtracted, B_singularity_subtraction_fixed
export B_finite_thickness
export get_curve
export plot_force_for_HSX, plot_integrand
export analytic_force_per_unit_length, interpolated_force_per_unit_length, force_finite_thickness, force_locally_circular_approximation
export force_finite_thickness_5D, force_finite_thickness_singularity_subtraction
export hifi_circular_coil_compute_Bz, hifi_circular_coil_force

include("utils.jl")
include("Curve.jl")
include("CurveCircle.jl")
include("CurveXYZFourier.jl")
include("configs.jl")
include("BiotSavart.jl")
include("results.jl")
include("force.jl")
include("hifi_circular_coil.jl")

end
