@doc """
	dyadic_rationals(I::DaubSupport, R::Int) -> Vector

The dyadic rationals of resolution `R` in the interval `I`.
"""->
function dyadic_rationals(I::DaubSupport, res::Int)
	@assert res >= 0
	collect( left(I):2.0^(-res):right(I) )
end

@doc """
	dyadic_rationals(I::DaubSupport, res, level) -> Vector

In a vector of dyadic rationals up to resolution `res` in `I`,
return the indices of those at exactly `level`.
"""->
function dyadic_rationals(I::DaubSupport, res::Int, level::Int)
	@assert 0 <= level <= res

	N = right(I)*2^res + 1

	if level == 0
		step = 2^res
		return 1:step:N
	else
		step = 2^(res-level)
		return 1+step:2*step:N
	end
end

@doc """
	dyadic_parent(i::Int, k::Int, L::Int)

In the vector `x` where the `i`'th entry is `x_i = (i-1)/2^L`,
`dyadic_parent` returns the index of `2x_i - k`.
"""->
function dyadic_parent(i::Int, k::Int, L::Int)
	2*i - 1 - k*2^L
end

@doc """
	isuniform(x::AbstractVector) -> Bool

Test if the elements in `x` increments with the same amount (with sign).
"""->
function isuniform(x::AbstractVector)
	@assert (Nx = length(x)) > 1

	if Nx == 2
		return true
	end

	diff = x[1] - x[2]
	for n = 3:Nx
		if !isapprox(diff, x[n-1] - x[n])
			return false
		end
	end

	return true
end

@doc """
	isdyrat(x::AbstractVector) -> Bool

Test if `x` is a vector of dyadic rationals.
"""->
function isdyrat(x::AbstractVector)
	res = -log2(x[2] - x[1])
	return isuniform(x) && isinteger( res )
end

