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

