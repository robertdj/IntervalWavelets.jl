using IntervalWavelets
using Test

function expected_boundary_function_values(side, R, index)
    value = if side == IntervalWavelets.LEFT
        if index == 0
            [0.0 ; 0.412364 ; -1.1266 ; 0.206182 ; 0.0] 
        elseif index == 1
            [0.0 ; 0.849475 ; 1.40586 ; -0.00765026 ; -0.365497 ; 0.0668906 ; 0.0]
        end
    elseif side == IntervalWavelets.RIGHT
        if index == 0
            [0.0 ; 0.445077 ; 0.651639 ; 0.890155 ; 0.0]
        elseif index == 1
            [0.0 ; 0.856098 ; 1.25341 ; 0.327042 ; 0.142971 ; -0.14055 ; 0.0]
        end
    end

    if R == 0
        return value[1:2:end]
    elseif R == 1
        return value
    end
end


@testset "Boundary scaling functions" begin
    @testset "Specific boundary scaling functions at $side and resolution $R" for side in [IntervalWavelets.LEFT, IntervalWavelets.RIGHT], R in 0:1
        h = interior_filter(2)
        phi = interior_scaling_function(h, 1)

        l = boundary_filters(2, side)

        phi = boundary_scaling_functions(l, phi, R)

        @test phi[0].(support(phi[0])) ≈ expected_boundary_function_values(side, R, 0) atol = 10.0^-4
        @test phi[1].(support(phi[1])) ≈ expected_boundary_function_values(side, R, 1) atol = 10.0^-4
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

