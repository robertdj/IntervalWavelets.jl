# ------------------------------------------------------------
# Dyadic rationals

"""
The value at the `k`'th index of Dyadic Rationals Vector of resolution 
`R` is the value of the vector/function at k/2^R.
"""
struct DyadicRationalsVector
	resolution::Int64
	values::OffsetVector{Float64, Vector{Float64}}

	# TODO: Check if resolution match length of values
	function DyadicRationalsVector(res, values) 
		if res < 0
			throw(DomainError("Resolution must be non-negative"))
		else
			new(res, values)
		end
	end
end

function DyadicRationalsVector(res, lower, upper)
	DyadicRationalsVector(
		res,
		OffsetVector(zeros(upper-lower+1), lower:upper)
		)
end

values(y::DyadicRationalsVector) = y.values

lower(y::DyadicRationalsVector) = linearindices(values(y))[1]
upper(y::DyadicRationalsVector) = linearindices(values(y))[end]
function is_in_support(k::Int64, y::DyadicRationalsVector)
	lower(y) <= k <= upper(y)
end

# TODO: checkindex instead of checkbounds?
function Base.checkbounds(y::DyadicRationalsVector, idx)
	checkbounds(y.values, idx)
end

function Base.getindex(y::DyadicRationalsVector, idx::Int)
	checkbounds(y, idx)
	y.values[idx]
end

function (y::DyadicRationalsVector)(k::Int64)
	if is_in_support(k, y)
		return values(y)[k]
	else
		return 0
	end
end

function index(y::DyadicRationalsVector)
	#= linspace(lower(y), upper(y), 2^y.res) =#
	lower(y):2.0^(-y.resolution):upper(y)
end

function Base.collect(y::DyadicRationalsVector)
	#= return x, y =#
end

function (y::DyadicRationalsVector)(k::Int64, res::Int64)
	if res < 0 
		error("Resolution must be non-negative")
	end
		
	if res > y.resolution
		error("Resolution is too large")
	end

	y(k*2^(y.resolution - res))
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

