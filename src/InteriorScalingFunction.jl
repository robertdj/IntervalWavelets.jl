struct InteriorScalingFunction
    values::OffsetArrays.OffsetVector{Float64}
    vanishing_moments::Int64
    scale::Int64
end


Base.values(phi::InteriorScalingFunction) = phi.values
vanishing_moments(phi::InteriorScalingFunction) = phi.vanishing_moments
scale(phi::InteriorScalingFunction) = phi.scale


"""
    interior_scaling_function(::InteriorFilter) -> InteriorScalingFunction

Compute an interior scaling function in the integers of its support.
"""
function interior_scaling_function(h::InteriorFilter)
    if vanishing_moments(h) == 1
        # Haar
        InteriorScalingFunction(OffsetArrays.OffsetVector([1.0; 0.0], 0:1), 1, 0)
    end

    H = dyadic_dilation_matrix(h)

    nonzero_function_values = eigval1(H)
    # Normalize scaling function in L2
    LinearAlgebra.rmul!(nonzero_function_values, 1/sum(nonzero_function_values))
    function_values = [0.0 ; nonzero_function_values ; 0.0]

    # TODO: We may need support and support boundaries
    support_begin, support_end = support(h)
    InteriorScalingFunction(
        OffsetArrays.OffsetVector(function_values, support_begin:support_end), 
        vanishing_moments(h), 0
    )
end


"""
    dyadic_dilation_matrix(::InteriorFilter) -> Matrix

Return the "dyadic dilation matrix" from an interior filter used to compute the associated scaling
function in the integers of its support.
"""
function dyadic_dilation_matrix(h::InteriorFilter)
    h_length = length(h)
    sz = h_length - 2

    H = Matrix{Float64}(undef, sz, sz)
    filter_coefficients = coefficients(h)

    left_support = support(h)[1]
    for col in 1:sz, row in 1:sz
        h_idx = 2*col - row + left_support
        H[col, row] = sqrt2 * h[h_idx]
    end

    return H
end

