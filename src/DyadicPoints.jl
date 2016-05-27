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
	Nx == 2 && return true

	const diff = x[1] - x[2]
	for n in 3:Nx
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
	@assert issorted(x)
	res = -log2(x[2] - x[1])
	res <= precision( eltype(x) ) && return false
	return isuniform(x) && isinteger( res )
end

#=
@doc """
	dyadic_rationals(dya_rat::Vector, level::Int)

In a vector of dyadic rationals return the indices of those up to `level`.
"""->
function dyadic_rationals(dy_rat::Vector{Float64}, level::Int)
	# Check input
	res = -log2(dy_rat[2])
	@assert isuniform(dy_rat) && isinteger( res ) "Input is not a vector of dyadic rationals"
	@assert 0 <= level <= res "Resolution must be greater than level"

	power2 = 2^level
	max_supp = dy_rat[end]
	Nlevel = Int(max_supp*power2+1)
	dyadic_level = Array{Int}(Nlevel)

	# An entry in dy_rat (at level res) is also in dyadic_level if
	# and only if it is an integer when multiplied with power2
	count = 1
	Nres = length(dy_rat)
	for n = 1:Nres
		if isinteger( power2*dy_rat[n] )
			dyadic_level[count] = n
			count += 1
		end
	end

	return dyadic_level
end

@doc """
	support(phi, J, k) -> lower, upper

From a vector of scaling function values `phi`, return the `lower` and
`upper` bound of the support of the version that is dilated with `J` and
translated with `k`.
"""->
function support(phi::Vector{Float64}, J::Int, k::Int)
	upper, res = factor_support( length(phi) )

	lower = k*2^(res-J) + 1
	upper = (upper+k)*2^(res-J) + 1

	return lower, upper
end

=#

#= function dyadic_level(y) =#
function dyadic_level(x)
	y = copy(x)
	level = 0
	while !isinteger(y)
		y *= 2
		level += 1
	end
	return level
end

