@doc """
	dyadic_rationals(I::DaubSupport, R::Int) -> Vector

The dyadic rationals of resolution `R` in the interval `I`.
"""->
function dyadic_rationals(I::DaubSupport, res::Integer)
	@assert res >= 0
	return left(I):2.0^(-res):right(I)
end

@doc """
	dyadic_rationals(I::DaubSupport, res, level) -> Vector

In a vector of dyadic rationals up to resolution `res` in `I`,
return the indices of those at exactly `level`.
"""->
function dyadic_rationals(I::DaubSupport, res::Integer, level::Integer)
	@assert 0 <= level <= res

	N = length(I)*2^res + 1

	if level == 0
		step = 2^res
		return 1:step:N
	else
		step = 2^(res-level)
		return 1+step:2*step:N
	end
end

@doc """
	isuniform(x::AbstractVector) -> Bool

Test if the elements in `x` increments with the same amount (with sign).
"""->
function isuniform(x::AbstractVector)
	@assert (Nx = length(x)) > 1
	Nx == 2 && return true

	diff = x[1] - x[2]
	for n in 3:Nx
		if !isapprox(diff, x[n-1] - x[n])
			return false
		end
	end

	return true
end

