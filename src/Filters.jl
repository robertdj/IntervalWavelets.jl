# --------------------------------------------------------------------------------------------------
# Functions and types for general filters

struct Filter
	coefficients::OffsetArrays.OffsetVector{Float64, Vector{Float64}}
end


function support(h::Filter)
    # TODO: Make this more high-level
    integers_in_support = h |> coefficients |> axes
    integers_in_support[1].indices
end


function support_boundaries(h::Filter)
    filter_support = support(h)

    return filter_support[1], filter_support[end]
end


coefficients(h::Filter) = h.coefficients
Base.length(h::Filter) = h |> coefficients |> length


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
		if p >= 0
			new(p, filter)
		else
			throw(DomainError(p, "Not a valid number of vanishing moments"))
		end
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

    return InteriorFilter(p, Filter(interior_filter(p, Val{phase})))
end


function interior_filter(p::Integer, ::Type{Val{:symmlet}})
    if !(2 <= p <= 8)
        throw(DomainError(p, "symmlet filter should have between 2 and 8 vanishing moments"))
    end

    OffsetArrays.OffsetVector(INTERIOR_COEFFICIENTS[p], -p+1:p)
end


function interior_filter(p::Integer, ::Type{Val{:min}})
    OffsetArrays.OffsetVector(Wavelets.wavelet(Wavelets.WT.Daubechies{p}()).qmf, 0:2*p-1)
end

