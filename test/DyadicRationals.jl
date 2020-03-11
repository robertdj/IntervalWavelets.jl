using IntervalWavelets
using Test


@testset "Dyadic Rationals" begin
    @testset "Construct Dyadic Rationals" begin 
        dr = DyadicRational(1, 2)

        @test IntervalWavelets.numerator(dr) == 1
        @test IntervalWavelets.scale(dr) == 2

        @test DyadicRational(4, 1) == DyadicRational(2, 0)
        @test DyadicRational(-4, 1) == DyadicRational(-2, 0)
        @test DyadicRational(5, 1) == DyadicRational(5, 1)

        @test_throws DomainError DyadicRational(1, -2)
    end


    @testset "Arithmetic with Dyadic Rationals" begin
        @test DyadicRational(2, 0) + 1 == DyadicRational(3, 0)
        @test DyadicRational(2, 1) + 1 == DyadicRational(4, 1)

        @test 2 * DyadicRational(2, 0) == DyadicRational(4, 0)
        @test 3 * DyadicRational(1, 2) == DyadicRational(3, 2)
    end


    @testset "All Dyadic Rationals" begin
        @test IntervalWavelets.all_dyadic_rationals(0, 1, 0) == DyadicRational.([0 ; 1], 0)
        @test IntervalWavelets.all_dyadic_rationals(0, 1, 1) == DyadicRational.(0:2, 1)
    end


    @testset "Convert Dyadic Rational" begin
        @test Integer(DyadicRational(3, 0)) == 3
        @test_throws InexactError Integer(DyadicRational(3, 2))

        @test float(DyadicRational(3, 0)) == 3.0
        @test float(DyadicRational(3, 2)) == 0.75
    end
end

