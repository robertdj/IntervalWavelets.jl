struct InteriorScalingFunction
    values::OffsetArrays.OffsetVector{Float64}
    filter::InteriorFilter
    # TODO: vanishing_moments is not needed when we have filter
    vanishing_moments::Int64
    scale::Int64
end


function InteriorScalingFunction(values, support, vanishing_moments, scale)
    InteriorScalingFunction(
        OffsetArrays.OffsetVector{Float64}(values, support[1]*2^scale:support[2]*2^scale),
        vanishing_moments, R
    )
end


Base.values(phi::InteriorScalingFunction) = phi.values
filter(phi::InteriorScalingFunction) = phi.filter
vanishing_moments(phi::InteriorScalingFunction) = phi.vanishing_moments
scale(phi::InteriorScalingFunction) = phi.scale


function support(phi::InteriorScalingFunction)
    #= DyadicRationalVector(axes(values(phi))[1].indices, scale(phi)) =#
    DyadicRational.(axes(values(phi))[1].indices, scale(phi))
end


# TODO: Does it make more sense to have phi(dr) instead of phi[dr]?
function Base.getindex(phi::InteriorScalingFunction, key::DyadicRational)
    key_scale = scale(key)
    phi_scale = scale(phi)

    if (phi_scale < key_scale)
        key = reduce(key)
        key_scale = scale(key)
        if (phi_scale < key_scale)
            throw(DomainError(key, "Scaling function not computed at this scale"))
        end
    end

    idx = numerator(key) * 2^(phi_scale - key_scale)
    if checkbounds(Bool, values(phi), idx)
        return values(phi)[numerator(key) * 2^(phi_scale - key_scale)]
    else
        return 0.0
    end
end


function Base.setindex!(phi::InteriorScalingFunction, value::Float64, key::DyadicRational)
    key_scale = scale(key)
    phi_scale = scale(phi)

    if (phi_scale < key_scale)
        key = reduce(key)
        key_scale = scale(key)
        if (phi_scale < key_scale)
            throw(DomainError(key, "Scaling function not computed at this scale"))
        end
    end
    #= if (phi_scale < key_scale) =#
    #=     throw(DomainError(key, "Scaling function not computed at this scale")) =#
    #= end =#

    idx = numerator(key) * 2^(phi_scale - key_scale)
    if checkbounds(Bool, values(phi), idx)
        values(phi)[idx] = value
    else
        throw(BoundsError())
    end
end


"""
    interior_scaling_function(::InteriorFilter) -> InteriorScalingFunction

Compute an interior scaling function in the integers in its support.
"""
function interior_scaling_function(h::InteriorFilter)
    # TODO: can InteriorFilter have 1 vanishing moment?
    if vanishing_moments(h) == 1
        # Haar
        InteriorScalingFunction(OffsetArrays.OffsetVector([1.0; 0.0], 0:1), 1, 0)
    end

    H = dyadic_dilation_matrix(h)

    nonzero_function_values = eigval1(H)
    # Normalize scaling function in L2
    LinearAlgebra.rmul!(nonzero_function_values, 1/sum(nonzero_function_values))
    function_values = [0.0 ; nonzero_function_values ; 0.0]

    InteriorScalingFunction(
        OffsetArrays.OffsetVector(function_values, support(h)), 
        h, vanishing_moments(h), 0
    )
end


"""
    interior_scaling_function(::InteriorFilter, ::Integer) -> InteriorScalingFunction

Compute an interior scaling function in the dyadic rationals of scale `R` in its support.
"""
function interior_scaling_function(h::InteriorFilter, R::Integer)
    if R < 0
        throw(DomainError(R, "Scale must be positive"))
    end

    phi = interior_scaling_function(h)

    for _ = 1:R
        phi = increase_resolution(phi)
    end

    return phi
end


"""
    increase_resolution(::InteriorScalingFunction)

Increase the resolution of a DaubScaling scaling function by one.
"""
function increase_resolution(phi::InteriorScalingFunction)
    # TODO: Function to get all dyadic rationals in interval
    support_left, support_right = support_boundaries(phi.filter)
    support2 = all_dyadic_rationals(support_left, support_right, scale(phi) + 1)

    R = scale(phi) + 1
    phi2 = InteriorScalingFunction(
        OffsetArrays.OffsetVector{Float64}(undef, support_left*2^R:support_right*2^R),
        phi.filter, vanishing_moments(phi), scale(phi) + 1
    )

    for (index, dr) in enumerate(support2)
        if isodd(index)
            phi2[dr] = phi[dr]
        end

        if iseven(index)
            phi_val = 0.0
            for j in support(phi.filter)
                phi_val += sqrt2 * phi.filter[j] * phi[2 * dr - j]
            end

            phi2[dr] = phi_val
        end
    end

    return phi2
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


function (phi::InteriorScalingFunction)(x::DyadicRational)
    phi[x]
end


function (phi::InteriorScalingFunction)(x::DyadicRational, J::Integer, k::Integer = 0)
end


@recipe function f(phi::InteriorScalingFunction)
    seriestype := :path

    x = phi |> support .|> float
    y = phi |> values |> parent

    @series begin
        x, y
    end
end

