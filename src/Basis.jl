struct IntervalScalingFunctionBasis
    left::BoundaryScalingFunctions{LeftScalingFunction}
    interior::InteriorScalingFunction
    right::BoundaryScalingFunctions{RightScalingFunction}

    scale::Int64

    left_boundary::DyadicRational
    right_boundary::DyadicRational

    function IntervalScalingFunctionBasis(left, interior, right, scale, left_boundary, right_boundary)
        if right_boundary <= left_boundary
            throw(DomainError(left_boundary, "Reconstruction interval is degenerate"))
        end

        p = vanishing_moments(left)
        if !(p == vanishing_moments(interior) == vanishing_moments(right))
            throw(DomainError(p, "All basis function must have the same number of vanishing moments"))
        end

        R = resolution(left)
        if R < scale
            throw(DomainError(scale, "Scale of basis functions must be greater than the resolution"))
        end

        length_of_support = 2*p - 1
        if length_of_support > 2^scale
            throw(DomainError(scale, "Number of vanishing moments is too large for the scale"))
        end

        if !(R == resolution(interior) == resolution(right))
            throw(DomainError(R, "All basis function must be evaluated at the same resolution"))
        end

        if R < resolution(left_boundary) || R < resolution(right_boundary)
            throw(DomainError(R, "The resolution of basis functions should be at least that of the endpoints"))
        end

        new(left, interior, right, scale, left_boundary, right_boundary)
    end
end


vanishing_moments(B::IntervalScalingFunctionBasis) = vanishing_moments(B.left)
resolution(B::IntervalScalingFunctionBasis) = resolution(B.left)
Base.size(B::IntervalScalingFunctionBasis) = 2^B.scale


# TODO: Should R have a default value based on J?
"""
    interval_scaling_function_basis(p, J, R; a, b)

Return a basis of scaling functions with `p` vanishing moments at scale `J` evaluated at resolution `R` on the interval `[a, b]`.
"""
function interval_scaling_function_basis(p::Integer, J::Integer, R::Integer; left_boundary = DyadicRational(0, 0), right_boundary = DyadicRational(1, 0))
    h = interior_filter(p)
    phi = interior_scaling_function(h, R)

    left = boundary_scaling_functions(LEFT, p, R)
    right = boundary_scaling_functions(RIGHT, p, R)

    IntervalScalingFunctionBasis(left, phi, right, J, left_boundary, right_boundary)
end


function Base.show(io::IO, B::IntervalScalingFunctionBasis)
    print(io, "Basis of scaling functions with $(vanishing_moments(B)) vanishing moments for [$(B.left_boundary), $(B.right_boundary)] at resolution $(resolution(B.left))")
end


function evaluate_function(B::IntervalScalingFunctionBasis, k::Integer)
    N = size(B)
    p = vanishing_moments(B)
    R = resolution(B)
    J = B.scale

    if 1 <= k <= p
        x_index = 0:2^(R - J)*(p + k)
        x = support(B.left[k - 1])
        y = B.left[k - 1].(x)
    elseif p < k <= N - p
        x_index = 2^(R - J)*(k - p + 1):2^(R - J)*(k + p)
        # TODO: This is the most widely used branch and we perform the same computation every time
        x = support(B.interior)
        y = B.interior.(x)
    elseif N - p < k <= N
        x_index = 2^(R - J)*(p + k):2^(R - J)*(p + k)
        x = support(B.right[N - k])
        y = B.right[N - k].(x)
    else 
        throw(DomainError(k, "There are only $N functions in this basis"))
    end

    return x_index, y
end


function reconstruct(coef::AbstractVector, B::IntervalScalingFunctionBasis)
    R = resolution(B)
    x = all_dyadic_rationals(0, 1, resolution(B))
    x_translated = x .- B.left_boundary

    N = length(x)
    #= y = Vector{Float64}(undef, N) =#
    y = zeros(Float64, N)

    p = vanishing_moments(B)

    s = supports(B.left)

    return x_translated, y
end

