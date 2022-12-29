using CoilForces
using Test

@testset "Regenerate the main results" begin
    @testset "Test plot_force_convergence_grid" begin
        CoilForces.plot_force_convergence_grid()
    end
end
