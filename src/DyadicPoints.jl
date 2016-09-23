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

typealias DaubSupport Interval{Int64}

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

@doc """
	dyadic_rationals(I::DaubSupport, R) -> Vector

The dyadic rationals of resolution `R` in the interval `I`.
"""->
function dyadic_rationals(I::DaubSupport, R::Integer)
	R >= 0 || throw(DomainError())
	return left(I):2.0^(-R):right(I)
end

@doc """
	dyadic_rationals(I::DaubSupport, R, level) -> Vector

In a vector of dyadic rationals up to resolution `R` in `I`,
return the **indices** of those at exactly `level`.
"""->
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

@doc """
	isuniform(x::AbstractVector) -> Bool

Test if the elements in `x` increments with the same amount (with sign).
"""->
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

