using IntervalWavelets
using Test


@testset "Scaling function basis" begin
    @testset "Non-default reconstruction intervals" begin
        J = 3
        R = 3
        l = DyadicRational(-1, 1)
        r = DyadicRational(1, 1)
        B = interval_scaling_function_basis(2, J, R; left_boundary = l, right_boundary = r)

        x = reconstruct(B, ones(2^J))[1]

        @test x == all_dyadic_rationals(l, r, R)
    end


    @testset "Reconstruct from unit vectors" begin
        p = 2
        J = 3
        R = 4
        B = interval_scaling_function_basis(p, J, R)

        x = all_dyadic_rationals(B.left_boundary, B.right_boundary, R)
        L = boundary_scaling_functions(LEFT, p, R)
        for k in 1:2
            coef = zeros(2^J)
            coef[k] = 1

            y = 2.0^(J/2) * L[k - 1].(2^J * x)
            xx, z = reconstruct(B, coef)

            @test x == xx
            @test y ≈ z atol = 1e-8
        end

        phi = interior_scaling_function(interior_filter(p), R)
        for k in 3:6
            coef = zeros(2^J)
            coef[k] = 1

            translation = k - 1
            y = 2.0^(J/2) * phi.(2^J * x .- translation)
            xx, z = reconstruct(B, coef)

            @test x == xx
            @test y ≈ z atol = 1e-8
        end

        right = boundary_scaling_functions(RIGHT, p, R)
        for k in 7:8
            coef = zeros(2^J)
            coef[k] = 1

            y = 2.0^(J/2) * right[8 - k].(2^J * (x .- B.right_boundary))
            xx, z = reconstruct(B, coef)

            @test x == xx
            @test y ≈ z atol = 1e-8
        end
    end


    @testset "Reconstruct constant function" begin
        B = interval_scaling_function_basis(2, 3, 4)

        # Coefficients found manually by solving [ B[j](x_i) ] * coef = ones()
        coef = [
            0.128004734,
            0.354064418,
            0.353553390,
            0.353553390,
            0.353553390,
            0.353553390,
            0.385317706,
            0.458021222
        ]

        y = reconstruct(B, coef)[2]

        @test y ≈ ones(17) atol = 1e-8
    end
end
