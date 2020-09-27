struct IntervalScalingFunctionBasis
    left::Vector{Vector{Float64}}
    interior::Vector{Float64}
    right::Vector{Vector{Float64}}

    vanishing_moments::Int64
    scale::Int64
    resolution::Int64

    left_boundary::DyadicRational
    right_boundary::DyadicRational
    support::Vector{DyadicRational}

    function IntervalScalingFunctionBasis(left, interior, right, p, scale, R, left_boundary, right_boundary)
        if right_boundary <= left_boundary
            throw(error("Reconstruction interval is degenerate"))
        end

        if R < scale
            throw(error("Scale of basis functions must be greater than the R"))
        end

        length_of_support = 2*p - 1
        if length_of_support > 2^scale
            throw(DomainError(scale, "Number of vanishing moments is too large for the scale"))
        end

        left_lengths = map(length, left)
        right_lengths = map(length, right)
        if !all(length(interior) .== left_lengths .== right_lengths)
            throw(error("All basis functions must have the same length"))
        end

        if R < resolution(left_boundary) || R < resolution(right_boundary)
            throw(DomainError(R, "The R of basis functions should be at least that of the endpoints"))
        end

        s = all_dyadic_rationals(left_boundary, right_boundary, R)

        new(left, interior, right, p, scale, R, left_boundary, right_boundary, s)
    end
end


vanishing_moments(B::IntervalScalingFunctionBasis) = B.vanishing_moments
resolution(B::IntervalScalingFunctionBasis) = B.resolution
Base.size(B::IntervalScalingFunctionBasis) = 2^B.scale
n_eval(B::IntervalScalingFunctionBasis) = length(B.interior)


"""
    interval_scaling_function_basis(p, J, R; a, b)

Return a basis of scaling functions with `p` vanishing moments at scale `J` evaluated at resolution `R` on the interval `[a, b]`.
"""
function interval_scaling_function_basis(p::Integer, J::Integer, R::Integer; left_boundary = DyadicRational(0, 0), right_boundary = DyadicRational(1, 0))
    h = interior_filter(p)
    phi = interior_scaling_function(h, R - J)
    interior_values = 2.0^(J/2) * phi.(support(phi))

    # Ensure that the boundary values are computed
    resolution = max(R - J, ceil(Int, log2(p)) + 1)

    left = boundary_scaling_functions(LEFT, p, resolution)
    left_x = all_dyadic_rationals(support_boundaries(left[p - 1])..., R - J)
    left_values = [2.0^(J/2) * left[k].(left_x) for k in 0:p - 1]

    right = boundary_scaling_functions(RIGHT, p, resolution)
    right_x = all_dyadic_rationals(support_boundaries(right[p - 1])..., R - J)
    right_values = [2.0^(J/2) * right[k].(right_x) for k in 0:p - 1]

    IntervalScalingFunctionBasis(left_values, interior_values, right_values, p, J, R, left_boundary, right_boundary)
end


function Base.show(io::IO, B::IntervalScalingFunctionBasis)
    print(io, "Basis of scaling functions with $(vanishing_moments(B)) vanishing moments for [$(B.left_boundary), $(B.right_boundary)] at resolution $(resolution(B))")
end


function Base.getindex(B::IntervalScalingFunctionBasis, k::Integer)
    N = size(B)
    p = vanishing_moments(B)
    R = resolution(B)
    J = B.scale

    y = Vector{Float64}(undef, n_eval(B))
    x_index = evaluate!(y, B, k)

    if 1 <= k < p
        y = y[1:length(x_index)]
    elseif N - p + 1 < k <= N
        start_index = length(y) - length(x_index) + 1
        y = y[start_index:end]
    end

    return B.support[x_index], y
end


function evaluate!(y::Vector{Float64}, B::IntervalScalingFunctionBasis, k::Integer)
    N = size(B)
    p = vanishing_moments(B)
    R = resolution(B)
    J = B.scale

    if 1 <= k <= p
        x_index = 1:(2p - 1)*2^(R - J) + 1
        copy!(y, B.left[k])
    elseif p < k <= N - p
        x_index = (k - p)*2^(R - J) + 1:(p + k - 1)*2^(R - J) + 1
        copy!(y, B.interior)
    elseif N - p < k <= N
        x_index = (N - 2p + 1)*2^(R - J) + 1:2^R + 1
        copy!(y, B.right[N - k + 1])
    else 
        throw(DomainError(k, "Basis functions are indexed from 1 through $N"))
    end

    return x_index
end


function reconstruct(B::IntervalScalingFunctionBasis, coef::AbstractVector)
    if length(coef) != size(B)
        throw(DimensionMismatch("The number of coefficients does not match the number of basis functions"))
    end

    x = all_dyadic_rationals(0, 1, resolution(B))
    x_translated = (B.right_boundary - B.left_boundary) * x .+ B.left_boundary

    y = Vector{Float64}(undef, n_eval(B))
    reconstruction = zeros(Float64, length(x))

    for (index, value) in enumerate(coef)
        x_index = evaluate!(y, B, index)
        reconstruction[x_index] += value * y
    end

    return x_translated, reconstruction
end


function reconstruct(B::IntervalScalingFunctionBasis, coef::AbstractMatrix)
    if !all(size(coef) .== size(B))
        throw(DimensionMismatch("The number of coefficients does not match the number of basis functions"))
    end

    N = n_eval(B)
    row_values = Vector{Float64}(undef, n_eval(B))
    col_values = similar(row_values)
    y = col_values * row_values'

    x = all_dyadic_rationals(0, 1, resolution(B))
    reconstruction = zeros(Float64, length(x), length(x))

    for row_count in 1:size(coef, 1)
        row_index = evaluate!(row_values, B, row_count)
        for col_count in 1:size(coef, 2)
            col_index = evaluate!(col_values, B, col_count)

            fill!(y, 0.0)
            alpha = coef[col_count, row_count]
            LinearAlgebra.BLAS.ger!(alpha, col_values, row_values, y)

            reconstruction[col_index, row_index] += y
        end
    end

    return reconstruction
end


@recipe function f(B::IntervalScalingFunctionBasis)
    for idx in 1:size(B)
        @series begin
            x, y = B[idx]
            float.(x), y
        end
    end
end

