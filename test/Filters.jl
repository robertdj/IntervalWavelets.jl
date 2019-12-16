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

    @testset "Construct interior filters" begin
        for phase in [:min; :symmlet]
            h = interior_filter(p, phase)
            @test vanishing_moments(h) == p
            @test length(h) == 2*p
        end
    end

    @testset "Daubechies filters are normalized in l2" begin
        for phase in [:min; :symmlet]
            h = interior_filter(p, phase)
            @test sum(coefficients(h).^2) ≈ 1 atol = sqrt(eps())
        end
    end

    @testset "Failure constructing interior filters" begin
        @test_throws DomainError interior_filter(0, :symmlet)
        @test_throws DomainError interior_filter(9, :symmlet)

        @test_throws MethodError interior_filter(p, :foo)
    end
end

