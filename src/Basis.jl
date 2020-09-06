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


# TODO: Should R have a default value based on J?
function interval_scaling_function_basis(p::Integer, J::Integer, R::Integer; left_boundary = DyadicRational(0, 0), right_boundary = DyadicRational(1, 0))
    h = interior_filter(p)
    phi = interior_scaling_function(h, R)

    left = boundary_scaling_functions(LEFT, p, R)
    right = boundary_scaling_functions(RIGHT, p, R)

    IntervalScalingFunctionBasis(left, phi, right, J, left_boundary, right_boundary)
end


function Base.show(io::IO, B::IntervalScalingFunctionBasis)
    print(io, "Basis of scaling functions with ", vanishing_moments(B), " vanishing moments for [", 
               B.left_boundary, ", ", B.right_boundary, "] at resolution ", resolution(B.left))
end


function reconstruct(coef, B::IntervalScalingFunctionBasis)
end

