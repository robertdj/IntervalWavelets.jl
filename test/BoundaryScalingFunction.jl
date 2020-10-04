using IntervalWavelets
using Test

# Side, index, resolution
expected_boundary_function_values = Dict{Tuple{Sides, Int, Int}, Vector{Float64}}(
    (LEFT, 0, 0) => [0.0 ; -1.1266 ; 0.0],
    (LEFT, 1, 0) => [0.0 ; 1.40586 ; -0.365497 ; 0.0],
    (LEFT, 0, 1) => [0.0 ; 0.412364 ; -1.1266 ; 0.206182 ; 0.0],
    (LEFT, 1, 1) => [0.0 ; 0.849475 ; 1.40586 ; -0.00765026 ; -0.365497 ; 0.0668906 ; 0.0],
    (RIGHT, 0, 0) => [0.0 ; 0.651639 ; 0.0],
    (RIGHT, 1, 0) => [0.0 ; 1.25341 ; 0.142971 ; 0.0],
    (RIGHT, 0, 1) => [0.0 ; 0.445077 ; 0.651639 ; 0.890155 ; 0.0],
    (RIGHT, 1, 1) => [0.0 ; 0.856098 ; 1.25341 ; 0.327042 ; 0.142971 ; -0.14055 ; 0.0]
)


@testset "$side boundary scaling functions" for side in [IntervalWavelets.LEFT, IntervalWavelets.RIGHT]
    @testset "Specific boundary scaling functions at resolution $R" for R in 0:1, index in 0:1
        phi = boundary_scaling_functions(side, 2, R)

        boundary_function_values = phi[index].(support(phi[index]))
        @test boundary_function_values ≈ expected_boundary_function_values[(side, index, R)] atol = 10.0^-4
    end


    @testset "Function values at 0" begin
        phi = boundary_scaling_functions(side, 2, 2)

        x = DyadicRational(0,0)
        if side == LEFT
            @test phi[0](x) ≈ 1.95133 atol = 10.0^-4
            @test phi[1](x) ≈ 0.293094 atol = 10.0^-4
        else
            @test phi[0](x) ≈ 1.128671 atol = 10.0^-4
            @test phi[1](x) ≈ -0.42407 atol = 10.0^-4
        end
    end


    @testset "Number of vanishing moments must be sufficiently high" begin
        @test_throws DomainError boundary_scaling_functions(side, 0, 0)
        @test_throws DomainError boundary_scaling_functions(side, 1, 0)
    end


    @testset "Increase resolution of boundary scaling function" for p in 2:3, index in 0:p - 1
        phi = boundary_scaling_functions(side, p, 0)

        phi0_support = support(phi[index])
        expected_support = IntervalWavelets.all_dyadic_rationals(support_boundaries(phi[index])..., 0)
        @test phi0_support == expected_support

        phi1 = IntervalWavelets.increase_resolution(phi)
        expected_support = IntervalWavelets.all_dyadic_rationals(support_boundaries(phi1[index])..., 1)
        @test support(phi1[index]) == expected_support

        @test phi1[index].(phi0_support) == phi[index].(phi0_support)
    end
end

