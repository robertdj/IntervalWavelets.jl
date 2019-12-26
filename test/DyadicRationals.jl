using IntervalWavelets
using Test


@testset "Dyadic Rationals" begin
    numerator = rand(-9:9)
    scale = rand(1:9)


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


    @testset "Arithmetic with Dyadic Rationals" begin
        @test DyadicRational(2, 0) + 1 == DyadicRational(3, 0)
        @test DyadicRational(2, 1) + 1 == DyadicRational(4, 1)

        @test IntervalWavelets.reduce(DyadicRational(4, 1)) == DyadicRational(2, 0)

        @test 2 * DyadicRational(2, 0) == DyadicRational(4, 0)
        @test 3 * DyadicRational(1, 2) == DyadicRational(3, 2)
    end


    @testset "All Dyadic Rationals" begin
        @test IntervalWavelets.all_dyadic_rationals(0, 1, 0) == DyadicRational.([0 ; 1], 0)
        @test IntervalWavelets.all_dyadic_rationals(0, 1, 1) == DyadicRational.(0:2, 1)
    end


    @testset "Convert Dyadic Rational" begin
        @test Integer(DyadicRational(3, 0)) == 3

        @test_throws InexactError Integer(DyadicRational(3, 1))
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
        @test_throws DomainError DyadicRationalVector(numerator, -scale - 1)

        @test_throws InexactError DyadicRationalVector(rand(2), scale)
        @test_throws InexactError DyadicRationalVector(numerator, rand())
        @test_throws InexactError DyadicRationalVector(rand(2), rand())
    end
end

