using IntervalWavelets
using OffsetArrays
using Test


@testset "Filters" begin
    filter_coefficients = OffsetVector(rand(3), -1:1)
    f = Filter(filter_coefficients)
    
    @test coefficients(f) == filter_coefficients
    @test support(f) == -1:1
    @test support_boundaries(f) == (-1, 1)
    @test length(f) == 3
    
    @test f[-1] == filter_coefficients[-1]
    @test f[0] == filter_coefficients[0]
    @test f[1] == filter_coefficients[1]
    
    @test f[2] == 0.0
    @test f[-2] == 0.0
end


@testset "Interior filters" begin
    p = rand(2:8)

    @testset "Construct interior filters $phase" for phase in [:min, :symmlet]
        h = interior_filter(p, phase)
        @test vanishing_moments(h) == p
        @test length(h) == 2*p
    end


    @testset "Daubechies filters are normalized in l2 $phase" for phase in [:min, :symmlet]
        h = interior_filter(p, phase)
        @test sum(coefficients(h).^2) ≈ 1 atol = sqrt(eps())
    end


    @testset "Failure constructing interior filters" begin
        @test_throws DomainError interior_filter(0, :symmlet)
        @test_throws DomainError interior_filter(9, :symmlet)

        @test_throws MethodError interior_filter(p, :foo)
    end
end


@testset "Boundary filters" begin
    @testset "Construct boundary filter" for s in [IntervalWavelets.LEFT, IntervalWavelets.RIGHT]
        p = 2

        b = boundary_filters(p, s)
        @test side(b) == s
        @test vanishing_moments(b) == p
        @test length(filters(b)) == p
    end


    @testset "Failure constructing boundary filter" begin
        @test_throws DomainError boundary_filters(0, IntervalWavelets.LEFT)
        @test_throws DomainError boundary_filters(9, IntervalWavelets.LEFT)
    end
end

