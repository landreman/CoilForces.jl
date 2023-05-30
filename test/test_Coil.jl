using CoilForces
using Test

@testset "Test winding pack angle functions" begin
    @testset "WindingPackAngleZero should always return 0" begin
        wpa = WindingPackAngleZero()
        n = 10
        for j in 1:n
            ϕ = (j - 1) * 2π / n
            @test get_winding_pack_angle(wpa, ϕ) ≈ 0
        end
    end

    @testset "WindingPackAngleFourier should return values consistent with the supplied Fourier series" begin
        cos_coeffs = [0.3, 0.7, -0.1, 1.1]
        sin_coeffs = [0.0, -0.2, 0.6, -0.4]
        nfourier = length(cos_coeffs)
        wpa = WindingPackAngleFourier(cos_coeffs, sin_coeffs)
        nϕ = 10
        for jϕ in 1:nϕ
            ϕ = (jϕ - 1) * 2π / nϕ

            should_be = 0.0
            for j = 1:nfourier
                m = j - 1
                sin_mϕ = sin(m * ϕ)
                cos_mϕ = cos(m * ϕ)
                should_be += cos_coeffs[j] * cos_mϕ + sin_coeffs[j] * sin_mϕ
            end

            @test get_winding_pack_angle(wpa, ϕ) ≈ should_be
        end
    end

end