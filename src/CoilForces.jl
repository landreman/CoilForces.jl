module CoilForces

import LinearAlgebra
using HCubature
using Plots
using Plots.PlotMeasures
using Printf
using Dates
using Interpolations
using CSV
import SpecialFunctions

export μ0, dot, cross, norm, normsq
export CurveCircle, CurveXYZFourier, fit_circle
export γ, γ_and_derivative, γ_and_2_derivatives, γ_and_3_derivatives, position_and_tangent, curvature, Frenet_frame, Frenet_frame_without_torsion, curve_length
export get_κ1_κ2, compute_regularization
export CoilCircularXSection, CoilRectangularXSection
export FrameCircle, FrameCentroid, FrameRotated, get_frame
export d_B_d_ϕ, B_filament_fixed, B_filament_adaptive
export d_B_d_ϕ_singularity_subtracted, B_singularity_subtraction_fixed, B_partition_fixed
export B_finite_thickness
export get_curve
export plot_force_for_HSX, plot_integrand
export analytic_force_per_unit_length, interpolated_force_per_unit_length, force_finite_thickness, force_locally_circular_approximation
export force_filament_adaptive
export force_finite_thickness_5D, force_finite_thickness_singularity_subtraction
export hifi_circular_coil_compute_Bz, hifi_circular_coil_force
export analytic_inductance_for_circular_coil, inductance_filament_adaptive, inductance_finite_thickness

include("utils.jl")
include("Curve.jl")
include("CurveCircle.jl")
include("CurveXYZFourier.jl")
include("configs.jl")
include("Coil.jl")
include("BiotSavart.jl")
include("results.jl")
include("force.jl")
include("hifi_circular_coil.jl")
include("vector_potential.jl")
include("inductance.jl")

end
