# ------------------------------------------------------------
# Dyadic rationals

"""
The value at the `k`'th index of Dyadic Rationals Vector of resolution 
`R` is k/2^R.
"""
struct DyadicRationalsVector
	# TODO: should values be any abstract vector?
	resolution::Int64
	parent::OffsetVector{Float64, Vector{Float64}}

	# TODO: Check if resolution match length of values
	function DyadicRationalsVector(res, values) 
		if res < 0
			throw(DomainError("Resolution must be non-negative"))
		else
			new(res, values)
		end
	end
end

function Base.show(io::IO, y::DyadicRationalsVector)
	println(io, "Dyadic rationals vector of resolution ", resolution(y),
				" and with indices ", linearindices(y))
	show(io, parent(y))
end

Base.parent(y::DyadicRationalsVector) = y.parent
resolution(y::DyadicRationalsVector) = y.resolution

Base.linearindices(y::DyadicRationalsVector) = linearindices(parent(y))

@inline function Base.getindex(y::DyadicRationalsVector, idx)
	parent(y)[idx]
end

@inline function Base.setindex!(y::DyadicRationalsVector, val, idx)
	parent(y)[idx] = val
end

function Base.length(y::DyadicRationalsVector)
	y |> linearindices |> length
end

function support(y::DyadicRationalsVector)
	2.0^(-y.resolution) * linearindices(y)
end

function Base.collect(y::DyadicRationalsVector)
	x = support(y)
	yvals = y |> parent |> parent

	return x, yvals
end

# Make a `plot` function for DyadicRationalsVector using the Plots package
@recipe function f(y::DyadicRationalsVector)
	seriestype := :path

	collect(y)
end

@recipe function f(Y::Vector{DyadicRationalsVector})
	for y in Y
		@series begin
			collect(y)
		end
	end
end


# ------------------------------------------------------------
# Intervals 

type Interval{T<:Real}
	left::T
	right::T

	function Interval(left, right)
		if left < right
			new(left, right)
		else
			throw(ArgumentError("Not an interval"))
		end
	end
end

function Interval(left, right)
	left, right = promote(left, right)
	T = typeof(left)
	Interval{T}(left, right)
end

const DaubSupport = Interval{Int64}

@inline left(I::Interval) = I.left
@inline right(I::Interval) = I.right
Base.length(I::Interval) = right(I) - left(I)
@compat Base.:+(S::DaubSupport, k::Integer) = DaubSupport(left(S)+k, right(S)+k)


function Base.intersect(I1::Interval, I2::Interval)
	if left(I1) >= right(I2) || left(I2) >= right(I1)
		return Void
	end

	L = max( left(I1), left(I2) )
	R = min( right(I1), right(I2) )

	DaubSupport(L, R)
end


# ------------------------------------------------------------
# Vectors on DaubSupport

"""
	dyadic_rationals(I::DaubSupport, R) -> Vector

The dyadic rationals of resolution `R` in the interval `I`.
"""
function dyadic_rationals(I::DaubSupport, R::Integer)
	R >= 0 || throw(DomainError())
	return left(I):2.0^(-R):right(I)
end

"""
	dyadic_rationals(I::DaubSupport, R, level) -> Vector

In a vector of dyadic rationals up to resolution `R` in `I`,
return the **indices** of those at exactly `level`.
"""
function dyadic_rationals(I::DaubSupport, R::Integer, level::Integer)
	0 <= level <= R || throw(DomainError())

	N = length(I)*2^R + 1

	if level == 0
		step = 2^R
		return 1:step:N
	else
		step = 2^(R-level)
		return 1+step:2*step:N
	end
end

"""
	isuniform(x::AbstractVector) -> Bool

Test if the elements in `x` increments with the same amount (with sign).
"""
function isuniform(x::AbstractVector)
	(Nx = length(x)) > 1 || throw(DomainError())
	Nx == 2 && return true

	diff = x[1] - x[2]
	@inbounds for n in 3:Nx
		if !isapprox(diff, x[n-1] - x[n])
			return false
		end
	end

	return true
end

