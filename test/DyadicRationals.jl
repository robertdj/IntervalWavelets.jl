using IntervalWavelets
using Test


@testset "Dyadic Rationals" begin
    numerator = rand(-9:9)
    scale = rand(1:9)

    @testset "Working Dyadic Rationals" begin 
        dr = DyadicRational(numerator, scale)

        @test IntervalWavelets.numerator(dr) == numerator
        @test IntervalWavelets.scale(dr) == scale
    end

    @testset "Wrong Dyadic Rationals" begin
        @test_throws DomainError DyadicRational(numerator, -scale)

        @test_throws InexactError DyadicRational(rand(), scale)
        @test_throws InexactError DyadicRational(numerator, rand())
        @test_throws InexactError DyadicRational(rand(), rand())
    end
end

