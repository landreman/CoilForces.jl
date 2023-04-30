using CoilForces
using Test

@testset "Regenerate the main results" begin
    @testset "Test plot_force_convergence_grid" begin
        CoilForces.plot_force_convergence_grid()
    end

    @testset "Test plot_modB_non_convergence_skipping_point_circular_with_ours" begin
        CoilForces.plot_modB_non_convergence_skipping_point_circular_with_ours()
    end

    @testset "Test plot_force_non_convergence_skipping_point_circular_with_ours" begin
        CoilForces.plot_force_non_convergence_skipping_point_circular_with_ours()
    end

    @testset "Test reproduce_Sienas_plot_of_locally_circular_approx" begin
        CoilForces.reproduce_Sienas_plot_of_locally_circular_approx()
    end
end
