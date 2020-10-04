# --------------------------------------------------------------------------------------------------
# Functions and types for general filters

struct Filter
	coefficients::OffsetArrays.OffsetVector{Float64, Vector{Float64}}
    support::AbstractRange{Int64}

    function Filter(coefficients::AbstractVector{Float64}, support::AbstractRange{Int64})
        new(OffsetArrays.OffsetVector(coefficients, support), support)
    end
end


coefficients(h::Filter) = h.coefficients
support(h::Filter) = h.support
Base.length(h::Filter) = h |> coefficients |> length


function support_boundaries(h::Filter)
    filter_support = support(h)

    return filter_support[1], filter_support[end]
end


function Base.show(io::IO, h::Filter)
	show(io, coefficients(h))
end


# TODO: How to handle indices with more than one integer outside support?
function Base.getindex(h::Filter, idx::Integer)
	if checkbounds(Bool, coefficients(h), idx)
		return coefficients(h)[idx]
	else
		return 0.0
	end
end


# --------------------------------------------------------------------------------------------------
# Functions and types for interacting with interior filters

struct InteriorFilter
	vanishing_moments::Int64
	filter::Filter

	function InteriorFilter(p, filter)
        if p < 0
			throw(DomainError(p, "Not a valid number of vanishing moments"))
        end

        new(p, filter)
	end
end


filter(h::InteriorFilter) = h.filter
vanishing_moments(h::InteriorFilter) = h.vanishing_moments
coefficients(h::InteriorFilter) = h |> filter |> coefficients
support(h::InteriorFilter) = h |> filter |> support
support_boundaries(h::InteriorFilter) = h |> filter |> support_boundaries
Base.length(h::InteriorFilter) = h |> filter |> length

Base.getindex(h::InteriorFilter, idx) = filter(h)[idx]


function Base.show(io::IO, h::InteriorFilter)
    left, right = support_boundaries(h)
	println(io, "Filter for Daubechies ", vanishing_moments(h), 
			" scaling function on [", left, ", ", right, "]:")
    show(io, filter(h))
end


"""
	interior_filter(p::Int[, phase=:symmlet])

Filter for Daubechies scaling function with `p` vanishing moments and `phase` being `:symmlet` or 
`:min` (minimum phase).
"""
function interior_filter(p::Integer, phase::Symbol=:symmlet)
	if p < 1 
		throw(DomainError(p, "Interior filter must have at least 1 vanishing moment"))
	end

    return InteriorFilter(p, interior_filter(p, Val{phase}))
end


function interior_filter(p::Integer, ::Type{Val{:symmlet}})
    if !(2 <= p <= 8)
        throw(DomainError(p, "symmlet filter should have between 2 and 8 vanishing moments"))
    end

    Filter(INTERIOR_COEFFICIENTS[p], -p+1:p)
end


function interior_filter(p::Integer, ::Type{Val{:min}})
    coefficients = Wavelets.wavelet(Wavelets.WT.Daubechies{p}()).qmf
    Filter(coefficients, 0:2*p-1)
end


# --------------------------------------------------------------------------------------------------
# Functions and types for interacting with boundary filters

struct BoundaryFilters
    side::Sides
    vanishing_moments::Int64
    filters::Vector{Filter}

    function BoundaryFilters(side, p, filters)
        if !(2 <= p <= 8)
            throw(DomainError(p, "Vanishing moments should be between 2 and 8"))
        end

        new(side, p, filters)
    end
end


side(B::BoundaryFilters) = B.side
vanishing_moments(B::BoundaryFilters) = B.vanishing_moments
filters(B::BoundaryFilters) = B.filters


function Base.getindex(B::BoundaryFilters, idx::Integer)
    filters(B)[idx + 1]
end


function boundary_filters(p::Integer, side::Sides)
    if !(2 <= p <= 8)
        throw(DomainError(p, "Vanishing moments should be between 2 and 8"))
    end

    supports = [0:(p + 2*k) for k in 0:(p - 1)]

    if side == LEFT
        coefficients = LEFT_SCALING_COEFFICIENTS[p]
    else
        coefficients = RIGHT_SCALING_COEFFICIENTS[p]
    end

    filters = map(Filter, coefficients, supports)

    return BoundaryFilters(side, p, filters)
end

