using IntervalWavelets
using OffsetArrays
using Test


@testset "Interior scaling functions" begin
    @testset "Dyadic dilation matrix has the expected form" begin
        h1 = interior_filter(3, :min)
        H1 = IntervalWavelets.dyadic_dilation_matrix(h1)

        @test size(H1) == (length(h1) - 2, length(h1) - 2)

        @test H1 == sqrt(2) * 
        [ h1[1] h1[0]    0     0  ;
          h1[3] h1[2] h1[1] h1[0] ;
          h1[5] h1[4] h1[3] h1[2] ;
             0     0  h1[5] h1[4] ]


        h2 = interior_filter(3, :symmlet)
        H2 = IntervalWavelets.dyadic_dilation_matrix(h2)

        @test H2 == sqrt(2) * 
        [ h2[-1] h2[-2]      0      0  ;
          h2[1]  h2[0]   h2[-1] h2[-2] ;
          h2[3]  h2[2]   h2[1]  h2[0] ;
             0      0    h2[3]  h2[2] ]
    end


    @testset "Construct interior scaling function" begin
        h = interior_filter(2)
        phi = interior_scaling_function(h)

        @test collect(support(h)) == Integer.(support(phi))
        @test vanishing_moments(phi) == 2
        @test scale(phi) == 0


        phi = interior_scaling_function(h, 1)

        @test vanishing_moments(phi) == 2
        @test scale(phi) == 1


        @test_throws DomainError interior_scaling_function(h, -1)
    end


    @testset "Specific interior scaling functions" begin
        h = interior_filter(2)
        phi = interior_scaling_function(h, 1)

        @test phi.(support(phi)) ≈ [0.0 ; 0.93301 ; 1.36602 ; 0.0 ; -0.36602 ; 0.06698 ; 0.0] atol = 10.0^-4
    end


    @testset "Haar scaling functions" begin
        haar_filter = interior_filter(1, :min)

        haar = interior_scaling_function(haar_filter)
        @test haar.(support(haar)) == [1.0 ; 0.0]

        haar = interior_scaling_function(haar_filter, 1)
        @test haar.(support(haar)) ≈ [1.0 ; 1.0 ; 0.0]
    end


    @testset "Increase resolution of interior scaling function" begin
        h = interior_filter(2)
        phi0 = interior_scaling_function(h)

        phi0_support = support(phi0)
        @test phi0_support == IntervalWavelets.all_dyadic_rationals(-1, 2, 0)

        phi1 = IntervalWavelets.increase_resolution(phi0)
        @test support(phi1) == IntervalWavelets.all_dyadic_rationals(-1, 2, 1)

        @test phi1.(phi0_support) == phi0.(phi0_support)
    end
end

