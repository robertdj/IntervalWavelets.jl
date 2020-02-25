using IntervalWavelets
using Test


@testset "Boundary scaling functions" begin
    @testset "Specific left boundary scaling functions" begin
        h = interior_filter(2)
        phi = interior_scaling_function(h, 1)

        l = boundary_filters(2, IntervalWavelets.LEFT)

        phi_left = boundary_scaling_functions(l, phi)

        @test phi_left[0].(support(phi_left[0])) ≈ [0.0 ; -1.1266 ; 0.0] atol = 10.0^-4
        @test phi_left[1].(support(phi_left[1])) ≈ [0.0 ; 1.40586 ; -0.365497 ; 0.0] atol = 10.0^-4
    end


    @testset "Specific right boundary scaling functions" begin
        h = interior_filter(2)
        phi = interior_scaling_function(h, 1)

        r = boundary_filters(2, IntervalWavelets.RIGHT)

        phi_right = boundary_scaling_functions(r, phi)

        @test phi_right[0].(support(phi_right[0])) ≈ [0.0 ; 0.6516 ; 0.0] atol = 10.0^-4
        @test phi_right[1].(support(phi_right[1])) ≈ [0.0 ; 1.2534 ; 0.14297 ; 0.0] atol = 10.0^-4
    end


    @testset "Increase resolution of boundary scaling function" begin
        #= h = interior_filter(2) =#
        #= phi0 = interior_scaling_function(h) =#

        #= phi0_support = support(phi0) =#
        #= @test phi0_support == IntervalWavelets.all_dyadic_rationals(-1, 2, 0) =#

        #= phi1 = IntervalWavelets.increase_resolution(phi0) =#
        #= @test support(phi1) == IntervalWavelets.all_dyadic_rationals(-1, 2, 1) =#

        #= @test phi1.(phi0_support) == phi0.(phi0_support) =#
    end
end

