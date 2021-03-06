abstract type AbstractScalingFunction end


struct InteriorScalingFunction <: AbstractScalingFunction
    values::OffsetArrays.OffsetVector{Float64, Vector{Float64}}
    support::Vector{DyadicRational}
    filter::InteriorFilter
    resolution::Int64
end


function Base.show(io::IO, phi::InteriorScalingFunction)
    print(io, "Interior scaling function with $(vanishing_moments(phi)) vanishing moments at resolution $(resolution(phi))")
end


Base.values(phi::AbstractScalingFunction) = phi.values
filter(phi::InteriorScalingFunction) = phi.filter
vanishing_moments(phi::InteriorScalingFunction) = phi |> filter |> vanishing_moments
resolution(phi::AbstractScalingFunction) = phi.resolution
support_boundaries(phi::InteriorScalingFunction) = support_boundaries(filter(phi))
support(phi::AbstractScalingFunction) = phi.support


function find_index(phi::AbstractScalingFunction, x::DyadicRational)
    x_scale = resolution(x)
    phi_scale = resolution(phi)

    if (phi_scale < x_scale)
        throw(DomainError(x, "Scaling function not computed at this resolution"))
    end

    idx = numerator(x) << (phi_scale - x_scale)
end


function get_value(phi::AbstractScalingFunction, x::DyadicRational)
    idx = find_index(phi, x)

    v = values(phi)
    if checkbounds(Bool, v, idx)
        return v[idx]
    else
        return 0.0
    end
end


function set_value!(phi::AbstractScalingFunction, x::DyadicRational, value::Float64)
    idx = find_index(phi, x)

    if checkbounds(Bool, values(phi), idx)
        values(phi)[idx] = value
    else
        throw(BoundsError())
    end
end


function set_value!(phi::AbstractScalingFunction, x::Integer, value::Float64)
    set_value!(phi, DyadicRational(x, 0), value)
end


function get_value(phi::AbstractScalingFunction, x::Integer)
    get_value(phi, DyadicRational(x, 0))
end


function (phi::InteriorScalingFunction)(x::DyadicRational)
    get_value(phi, x)
end


function (phi::InteriorScalingFunction)(x::DyadicRational, J::Integer, k::Integer = 0)
end


function initialize_interior_scaling_function(f::InteriorFilter, R::Integer)
    support_left, support_right = support_boundaries(f)
    support = all_dyadic_rationals(support_left, support_right, R)

    phi2 = InteriorScalingFunction(
        OffsetArrays.OffsetVector{Float64}(undef, support_left*2^R:support_right*2^R),
        support, f, R
    )
end


"""
    interior_scaling_function(::InteriorFilter) -> InteriorScalingFunction

Compute an interior scaling function in the integers in its support.
"""
function interior_scaling_function(h::InteriorFilter)
    if vanishing_moments(h) == 1
        # Haar
        return InteriorScalingFunction(
            OffsetArrays.OffsetVector([1.0; 0.0], 0:1), 
            DyadicRational.(0:1, 0), h, 0
        )
    end

    H = dyadic_dilation_matrix(h)

    nonzero_function_values = eigval1(H)
    # Normalize scaling function in L2
    LinearAlgebra.rmul!(nonzero_function_values, 1/sum(nonzero_function_values))
    function_values = [0.0 ; nonzero_function_values ; 0.0]

    InteriorScalingFunction(
        OffsetArrays.OffsetVector(function_values, support(h)), 
        DyadicRational.(support(h), 0), h, 0
    )
end


"""
    interior_scaling_function(::InteriorFilter, ::Integer) -> InteriorScalingFunction

Compute an interior scaling function in the dyadic rationals of scale `R` in its support.
"""
function interior_scaling_function(h::InteriorFilter, R::Integer)
    if R < 0
        throw(DomainError(R, "Resolution must be non-negative"))
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
    h = filter(phi)
    phi2 = initialize_interior_scaling_function(h, resolution(phi) + 1)
    support2 = support(phi2)

    for (index, x) in enumerate(support2)
        if isodd(index)
            set_value!(phi2, x, phi(x))
            continue
        end

        phi_val = 0.0
        for j in support(h)
            phi_val += h[j] * phi(2*x - j)
        end

        set_value!(phi2, x, sqrt2 * phi_val)
    end

    return phi2
end


"""
    dyadic_dilation_matrix(::InteriorFilter) -> Matrix

Return the "dyadic dilation matrix" from an interior filter used to compute the associated scaling
function in the integers of its support.
"""
function dyadic_dilation_matrix(h::InteriorFilter)
    sz = length(h) - 2
    left_support = support_boundaries(h)[1]

    [sqrt2 * h[2*i - j + left_support] for i = 1:sz, j = 1:sz]
end


function Base.collect(phi::AbstractScalingFunction)
    x = phi |> support .|> float
    y = phi |> values |> parent

    return x, y
end


@recipe function f(phi::AbstractScalingFunction)
    seriestype := :path

    collect(phi)
end

