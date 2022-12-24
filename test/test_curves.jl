using CoilForces
using Test

@testset "Test CurveXYZFourier basics" begin
    @testset "Test CurveXYZFourier constructor" begin
        xc = [1.1, 0.2]
        xs = [0.0, 3.1]
        yc = [0.9, -0.1]
        ys = [0.0, -0.3]
        zc = [-0.2, 1.4]
        zs = [-1.3, 0.2]
        c = CurveXYZFourier(xc, xs, yc, ys, zc, zs)
        @test c.n == 2
    end
end

@testset "Test Frenet frame functions" begin
    @testset "Test curvature and torsion for CurveCircle" begin
        R0 = 3.7
        c = CurveCircle(R0)
        ϕ = -0.3
        dℓdϕ, κ, τ = curvature_torsion(c, ϕ)
        @test dℓdϕ ≈ R0
        @test κ ≈ 1 / R0
        @test τ ≈ 0.0
    end

    @testset "For CurveXYZFourier, compare the curvature and torsion to reference values from simsopt" begin
        # Reference values for this test are from the simsopt calculation in 20221223-01-CurveXYZFourier_julia_simsopt_benchmark
        xc = [1.1, 0.2, -0.3]
        xs = [0.0, 3.1, -0.2]
        yc = [0.9, -0.1, -0.1]
        ys = [0.0, -0.3, 0.3]
        zc = [-0.2, 1.4, 0.2]
        zs = [-1.3, 0.2, 0.1]
        c = CurveXYZFourier(xc, xs, yc, ys, zc, zs)

        nϕ = 10
        ϕ = (collect(1:nϕ) .- 1) * 2π / nϕ
        γ_python = [[ 1.                 0.7                1.4              ]
            [ 2.801021279410141  0.897177980325815  1.207089893087926]
            [ 4.235226647243955  0.841018620799196  0.319810221738215]
            [ 4.346733950410966  0.550150868298702 -0.662994412970133]
            [ 2.857837088178225  0.488347469423712 -1.248368994420958]
            [ 0.600000000000001  0.9               -1.4              ]
            [-1.166854082553171  1.411652530576288 -1.293271792078916]
            [-1.784930551535976  1.473455929451277 -0.925859969029699]
            [-1.426209652869008  1.058981379200804 -0.178169435238341]
            [-0.462824678285131  0.679215221924206  0.781764488911906]]
        κ_python = [0.313157238338248 0.241466260874787 0.597096666515087 1.249424020105404 0.165639568171942 0.067127788602126 0.291779560176924 3.051815139115102 0.250500810930593 0.325566616134798]
        τ_python = [-0.311006548883502 -0.170657800876417  0.141523831415142 0.237320800055858  0.184397780420979 -0.17394423025195 -0.460068842056389  0.312083478401266  0.53525579692768 -0.298916628057598]
        for j in 1:nϕ
            position = γ(c, ϕ[j])
            @test position ≈ γ_python[j, :]
            dℓdϕ, κ, τ = curvature_torsion(c, ϕ[j])
            @test κ ≈ κ_python[j]
            @test τ ≈ τ_python[j]
        end
    end
end
