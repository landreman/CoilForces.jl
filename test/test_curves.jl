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

    @testset "Ensure γ_and_3_derivatives, γ_and_2_derivatives, γ_and_derivative, and γ are consistent" begin
        function subtest(curve)
            for ϕ in [j * 0.1 for j in 1:20]
                data4 = γ_and_3_derivatives(curve, ϕ)
                data3 = γ_and_2_derivatives(curve, ϕ)
                data2 = γ_and_derivative(curve, ϕ)
                data1 = γ(curve, ϕ)

                @test data1 ≈ data2[:, 1]
                @test data1 ≈ data3[:, 1]
                @test data1 ≈ data4[:, 1]

                @test data2 ≈ data3[:, 1:2]
                @test data2 ≈ data4[:, 1:2]

                @test data3 ≈ data4[:, 1:3]
            end
        end

        for j in 1:6
            curve = get_curve("hsx", j)
            subtest(curve)
        end
        R0 = 3.7
        curve = CurveCircle(R0)
        subtest(curve)
    end
end

@testset "Test Frenet frame functions" begin
    @testset "Test curvature and torsion for CurveCircle" begin
        R0 = 3.7
        c = CurveCircle(R0)
        for ϕ in range(-3π, 3π, length=50)
            dℓdϕ, κ, τ, γ0, tangent, normal, binormal = Frenet_frame(c, ϕ)
            @test dℓdϕ ≈ R0
            @test κ ≈ 1 / R0
            @test τ ≈ 0.0
            @test tangent ≈ [-sin(ϕ), cos(ϕ), 0]
            @test normal ≈ [-cos(ϕ), -sin(ϕ), 0]
            @test binormal ≈ [0, 0, 1]

            dℓdϕ, κ, γ0, tangent, normal = Frenet_frame_without_torsion(c, ϕ)
            @test dℓdϕ ≈ R0
            @test κ ≈ 1 / R0
            @test tangent ≈ [-sin(ϕ), cos(ϕ), 0]
            @test normal ≈ [-cos(ϕ), -sin(ϕ), 0]

            dℓdϕ, κ = curvature(c, ϕ)
            @test dℓdϕ ≈ R0
            @test κ ≈ 1 / R0

            position, tangent = position_and_tangent(c, ϕ)
            @test position ≈ [R0 * cos(ϕ), R0 * sin(ϕ), 0]
            @test tangent ≈ [-sin(ϕ), cos(ϕ), 0]
        end
    end

    @testset "For CurveXYZFourier, compare the curvature and torsion to reference values from simsopt" begin
        # Reference values for this test are from the simsopt calculation in 
        # 20221223-01-CurveXYZFourier_julia_simsopt_benchmark
        xc = [1.1, 0.2, -0.3, 0]
        xs = [0.0, 3.1, -0.2, 0]
        yc = [0.9, -0.1, -0.1, 0]
        ys = [0.0, -0.3, 0.3, 0]
        zc = [-0.2, 1.4, 0.2, 0]
        zs = [-1.3, 0.2, 0.1, 0]
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
        
        tangent_python = [[ 9.832820049844601e-01  1.092535561093844e-01  1.456714081458459e-01]
            [ 9.433185873360874e-01  6.373015368373459e-02 -3.257123121678022e-01]
            [ 6.460254166836512e-01 -1.634961935243430e-01 -7.456005335980960e-01]
            [-6.479891618322098e-01 -2.285171887199828e-01 -7.265603489095080e-01]
            [-9.830939584493932e-01  8.785804691431921e-02 -1.606462961063964e-01]
            [-9.684929325968498e-01  2.490410398106185e-01 -2.178437749328550e-17]
            [-9.474701204903810e-01  2.728021343545771e-01  1.669711540043080e-01]
            [-9.602967783940103e-02 -3.887383559867874e-01  9.163300669293519e-01]
            [ 5.575594309019976e-01 -3.939825528547506e-01  7.306881886675699e-01]
            [ 7.978178480950131e-01 -1.265279052827657e-01  5.894721116777266e-01]]

        normal_python = [[ 0.12475082096959   0.178560495210413 -0.975988412952544]
            [-0.286139414889231 -0.341046642167714 -0.89543923474068 ]
            [-0.744917181383399 -0.348207955763647 -0.569077861476541]
            [-0.760022677744787  0.256316443239754  0.597216384770683]
            [-0.070200770905801  0.629459546677168  0.773855626626314]
            [ 0.181579073184085  0.70614084016033   0.68439327439658 ]
            [ 0.0714347403468   -0.328365376623448  0.941845665332982]
            [ 0.949121578511439 -0.313136101108239 -0.033376809114799]
            [ 0.756410089792438  0.603752914041733 -0.251646964707129]
            [ 0.521309475169573  0.635926207462868 -0.569064398605568]]
        
        binormal_python = [[-0.13264136361314   0.977844471287569  0.16194585092302 ]
            [-0.16814957042198   0.9378836043868   -0.3034799278195  ]
            [-0.16658197345009   0.923049410511676 -0.346742313361841]
            [ 0.049755155126212  0.939192086521409 -0.339768522900414]
            [ 0.169109788671543  0.772050285078229 -0.612650174802955]
            [ 0.170442012695118  0.662830049369904 -0.729113054306894]
            [ 0.31176505360532   0.904298167048087  0.291628833320615]
            [ 0.299910870389881  0.866503475334877  0.399030320971432]
            [-0.342010609590457  0.693008056830633  0.634640509340047]
            [-0.302858238057786  0.761307131088793  0.573313474282513]]

        for j in 1:nϕ
            position = γ(c, ϕ[j])
            @test position ≈ γ_python[j, :]

            dℓdϕ, κ, τ, γ0, tangent, normal, binormal = Frenet_frame(c, ϕ[j])
            @test κ ≈ κ_python[j]
            @test τ ≈ τ_python[j]
            @test γ0 ≈ γ_python[j, :]
            @test tangent ≈ tangent_python[j, :]
            @test normal ≈ normal_python[j, :]
            @test binormal ≈ binormal_python[j, :] 

            dℓdϕ, κ, γ0, tangent, normal = Frenet_frame_without_torsion(c, ϕ[j])
            @test κ ≈ κ_python[j]
            @test γ0 ≈ γ_python[j, :]
            @test tangent ≈ tangent_python[j, :]
            @test normal ≈ normal_python[j, :]

            dℓdϕ2, κ2 = curvature(c, ϕ[j])
            @test κ2 ≈ κ_python[j]
            @test κ2 ≈ κ
            @test dℓdϕ2 ≈ dℓdϕ

            γ0, tangent = position_and_tangent(c, ϕ[j])
            @test γ0 ≈ γ_python[j, :]
            @test tangent ≈ tangent_python[j, :]
        end
    end
end

@testset "Test fit_circle" begin
    @testset "Fitting a circle to a circle, all aspects of the Frenet frame should be identical" begin
        curve1 = CurveCircle(1.7)
        for ϕ_fit in range(0, 10, length=7)
            curve2 = fit_circle(curve1, ϕ_fit)
            for ϕ_test in range(0, 10, length=8)
                differential_arclength1, curvature1, torsion1, position1, tangent1, normal1, binormal1 = Frenet_frame(curve1, ϕ_test)
                differential_arclength2, curvature2, torsion2, position2, tangent2, normal2, binormal2 = Frenet_frame(curve2, ϕ_test)
                @test differential_arclength1 ≈ differential_arclength2
                @test curvature1 ≈ curvature2
                @test torsion1 ≈ torsion2
                @test position1 ≈ position2
                @test tangent1 ≈ tangent2
                @test normal1 ≈ normal2
                @test binormal1 ≈ binormal2
            end
        end
    end

    @testset "Fitting a circle to a general curve, a few quantities should be identical at the fit point" begin
        curve1 = get_curve("hsx", 1)
        for ϕ in range(0, 10, length=7)
            curve2 = fit_circle(curve1, ϕ)
            differential_arclength1, curvature1, torsion1, position1, tangent1, normal1, binormal1 = Frenet_frame(curve1, ϕ)
            differential_arclength2, curvature2, torsion2, position2, tangent2, normal2, binormal2 = Frenet_frame(curve2, ϕ)
            @test curvature1 ≈ curvature2
            @test position1 ≈ position2
            @test tangent1 ≈ tangent2
            @test normal1 ≈ normal2
            @test binormal1 ≈ binormal2
        end
    end

end

@testset "Test curve utilities" begin
    @testset "Test curve_length for a circle" begin
        R = 17.2
        c = CurveCircle(R)
        @test curve_length(c) ≈ 2π * R
    end

    @testset "Test curve_length for an HSX coil" begin
        c = get_curve("hsx", 1)
        # Compare to a reference value from simsopt:
        @test curve_length(c) ≈ 2.054316451786527
    end

    @testset "Centroid of a CurveCircle should be the origin" begin
        c = CurveCircle(2.7)
        @test CoilForces.centroid(c) ≈ zeros(3) atol=1e-14
    end

    @testset "Centroid of an ellipse should equal the 0-frequency modes" begin
        x0 = 1.7
        y0 = -0.3
        z0 = -0.9

        curve = CurveXYZFourier([x0, 0.2], [0, -0.1], [y0, -1.4], [0, 0.3], [z0, 1.1], [0, -0.2])
        @test CoilForces.centroid(curve) ≈ [x0, y0, z0]
    end

    @testset "Check centroid of HSX coil against reference values from simsopt" begin
        curve = get_curve("hsx", 1)
        @test CoilForces.centroid(curve) ≈ [1.449921164520124, 0.080852101200521, 0.053270789528036]
    end
end