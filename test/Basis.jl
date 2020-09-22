using IntervalWavelets
using Test


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
    for k in 1:2
        coef = zeros(2^J)
        coef[k] = 1

        y = B.left[k - 1].(2^J * x)
        xx, z = reconstruct(B, coef)

        @test x == xx
        @test y == z
    end

    # for k in 3:6
    #     coef = zeros(2^J)
    #     coef[k] = 1

    #     y = B.interior.(2^J * x .- k)
    #     xx, z = reconstruct(B, coef)

    #     @test x == xx
    #     @test y == z
    # end

    for k in 7:8
        coef = zeros(2^J)
        coef[k] = 1

        y = B.right[8 - k].(2^J * x .- B.right_boundary)
        xx, z = reconstruct(B, coef)

        @test x == xx
        # @test y == z
    end
end

