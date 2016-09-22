type DyadicRational
	numerator::Int64
	denominator::Int64
	resolution::Int64
	float::Float64

	function DyadicRational(numerator, denominator, resolution, float)
		if resolution >= 0 && isequal(numerator/denominator, float) &&
			denominator == 2^resolution
			new(numerator, denominator, resolution, float)
		else
			throw(ArgumentError("Not a valid dyadic rational"))
		end
	end
end

DyadicRational(numerator::Integer, resolution::Integer) = DyadicRational(numerator, 2^resolution, resolution)
DyadicRational(numerator::Integer, denominator::Integer, resolution::Integer) = DyadicRational(numerator, denominator, resolution, numerator/denominator)

DyadicRational(x::Integer) = DyadicRational(x, 1, 0, float(x))
function DyadicRational(x::Real)
	num = x
	R = 0
	prec = precision(x)
	while !isinteger(num)
		if R >= prec
			throw(DomainError("Input is not a dyadic rational"))
		end

		R += 1
		num *= 2
	end

	DyadicRational(num, 2^R, R, x)
end

for name in [:+, :-]
	@eval begin
		Base.$name(x::DyadicRational, y::DyadicRational) = DyadicRational( $name(x.float, y.float) )
	end
end
Base.isless(x::DyadicRational, y::DyadicRational) = isless(x.float, y.float)

function Base.show(io::IO, x::DyadicRational)
	show(io, x.float)
end


@doc """
	dyadic_rationals(I::DaubSupport, R::Int) -> Vector

The dyadic rationals of resolution `R` in the interval `I`.
"""->
function dyadic_rationals(I::DaubSupport, R::Integer)
	R >= 0 || throw(DomainError())
	return left(I):2.0^(-R):right(I)
end

@doc """
	dyadic_rationals(I::DaubSupport, res, level) -> Vector

In a vector of dyadic rationals up to resolution `res` in `I`,
return the indices of those at exactly `level`.
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

