using IntervalWavelets
using Test


@testset "Dyadic Rationals" begin
    numerator = rand(-9:9)
    scale = rand(0:9)

    @testset "Construct Dyadic Rationals" begin 
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


@testset "Dyadic Rational vector" begin
    numerator = rand(-9:9, 2)
    scale = rand(0:9)

    @testset "Construct Dyadic Rational vector" begin
        drv = DyadicRationalVector(numerator, scale)

        @test IntervalWavelets.numerator(drv) == numerator
        @test IntervalWavelets.scale(drv) == scale
    end

    @testset "Use Dyadic Rational vector" begin
        drv = DyadicRationalVector(numerator, scale)

        @test isa(drv[1], DyadicRational)
        @test isa(drv[end], DyadicRational)
    end

    @testset "Wrong Dyadic Rationals vector" begin
        @test_throws DomainError DyadicRationalVector(numerator, -scale)

        @test_throws InexactError DyadicRationalVector(rand(2), scale)
        @test_throws InexactError DyadicRationalVector(numerator, rand())
        @test_throws InexactError DyadicRationalVector(rand(2), rand())
    end
end

