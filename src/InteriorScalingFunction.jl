abstract type AbstractScalingFunction end


struct InteriorScalingFunction <: AbstractScalingFunction
    values::OffsetArrays.OffsetVector{Float64, Vector{Float64}}
    filter::InteriorFilter
    # TODO: vanishing_moments is not needed when we have filter
    vanishing_moments::Int64
    scale::Int64
end


function InteriorScalingFunction(values, support, vanishing_moments, scale)
    InteriorScalingFunction(
        OffsetArrays.OffsetVector{Float64}(values, support[1]*2^scale:support[2]*2^scale),
        vanishing_moments, scale
    )
end


Base.values(phi::AbstractScalingFunction) = phi.values
filter(phi::InteriorScalingFunction) = phi.filter
vanishing_moments(phi::AbstractScalingFunction) = phi.vanishing_moments
scale(phi::AbstractScalingFunction) = phi.scale
support_boundaries(phi::InteriorScalingFunction) = support_boundaries(filter(phi))


function support(phi::AbstractScalingFunction)
    # TODO: Prettify
    #= DyadicRationalVector(axes(values(phi))[1].indices, scale(phi)) =#
    DyadicRational.(axes(values(phi))[1].indices, scale(phi))
end


function find_index(phi::AbstractScalingFunction, key::DyadicRational)
    key_scale = scale(key)
    phi_scale = scale(phi)

    if (phi_scale < key_scale)
        throw(DomainError(key, "Scaling function not computed at this scale"))
    end

    idx = numerator(key) << (phi_scale - key_scale)
end


function Base.getindex(phi::AbstractScalingFunction, key::DyadicRational)
    idx = find_index(phi, key)

    if checkbounds(Bool, values(phi), idx)
        return values(phi)[idx]
    else
        return 0.0
    end
end


function Base.setindex!(phi::AbstractScalingFunction, value::Float64, key::DyadicRational)
    idx = find_index(phi, key)

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
    if vanishing_moments(h) == 1
        # Haar
        return InteriorScalingFunction(OffsetArrays.OffsetVector([1.0; 0.0], 0:1), h, 1, 0)
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
    # TODO: Function for initializing phi?
    support_left, support_right = support_boundaries(filter(phi))
    R = scale(phi) + 1
    support2 = all_dyadic_rationals(support_left, support_right, R)

    phi2 = InteriorScalingFunction(
        OffsetArrays.OffsetVector{Float64}(undef, support_left*2^R:support_right*2^R),
        filter(phi), vanishing_moments(phi), R
    )

    h = filter(phi)

    for (index, dr) in enumerate(support2)
        if isodd(index)
            phi2[dr] = phi[dr]
        else
        #= elseif iseven(index) =#
            phi_val = 0.0
            for j in support(h)
                phi_val += sqrt2 * h[j] * phi[2*dr - j]
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

    left_support = support_boundaries(h)[1]
    for col in 1:sz, row in 1:sz
        h_idx = 2*col - row + left_support
        H[col, row] = sqrt2 * h[h_idx]
    end

    return H
end


# TODO: Does it make more sense to have only phi(dr) instead of phi[dr]?
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

