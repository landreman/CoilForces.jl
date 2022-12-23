using CoilForces
using Test

# Compare to reference values from 
# "20221019-02 Convergence of julia calculation of circular coil self-force.docx"
# and
# 20221016_04_Finite_thickness_circular_coil.jl

@testset "Test force for a circular coil" begin
    @testset "Test analytic solution" begin
        # Major radius of coil [meters]
        R0 = 2.3

        # Minor radius of coil [meters]
        a = 0.1

        # Total current [Amperes]
        I = 3.1e6

        curve = CurveCircle(R0)
        coil = Coil(curve, I, a)
        @test analytic_force_per_unit_length(coil) ≈ 1.8655666361140157e6
    end
end
