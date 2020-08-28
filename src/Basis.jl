struct IntervalScalingFunctionBasis
    left::BoundaryScalingFunctions{LeftScalingFunction}
    interior::InteriorScalingFunction
    right::BoundaryScalingFunctions{RightScalingFunction}

    left_boundary::DyadicRational
    right_boundary::DyadicRational

    function IntervalScalingFunctionBasis(left, interior, right, left_boundary, right_boundary)
        p = vanishing_moments(left)
        if !(p == vanishing_moments(interior) == vanishing_moments(right))
            throw(DomainError(p, "All basis function must have the same number of vanishing moments"))
        end

        R = resolution(left)
        if !(R == resolution(interior) == resolution(right))
            throw(DomainError(R, "All basis function must be evaluated at the same resolution"))
        end

        if R > resolution(left_boundary)
        end

        new(left, interior, right, left_boundary, right_boundary)
    end
end


vanishing_moments(B::IntervalScalingFunctionBasis) = vanishing_moments(B.left)


function interval_scaling_function_basis(p::Integer, R::Integer; left_boundary = DyadicRational(0, 0), right_boundary = DyadicRational(1, 0))
    h = interior_filter(p)
    phi = interior_scaling_function(h, R)

    left = boundary_scaling_functions(LEFT, p, R)
    right = boundary_scaling_functions(RIGHT, p, R)

    IntervalScalingFunctionBasis(left, phi, right, left_boundary, right_boundary)
end


function Base.show(io::IO, B::IntervalScalingFunctionBasis)
    print(io, "Basis of scaling functions with ", vanishing_moments(B), " vanishing moments for [", 
               B.left_boundary, ", ", B.right_boundary, "] at resolution ", resolution(B.left))
end


