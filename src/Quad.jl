@doc """
	trapezquad(x::Vector, y::Vector) -> Number

Approximate the integral of `y` over the the interval `[min(x), max(x)]`, using the divisions in `x`.
"""->
function trapezquad(x::AbstractVector, y::AbstractVector)
	@assert 2 <= (Nx = length(x)) == length(y)
	@assert issorted(x)

	integral = 0.0
	for n in 2:Nx
		@inbounds integral += (x[n]-x[n-1])*0.5*(y[n]+y[n-1])
	end

	return integral
end

@doc """
	inner(x::Vector, y::Vector, z::Vector) -> Number

Approximate the inner product between `y` and `z` over the the interval `[min(x), max(x)]`, using the divisions in `x`.
The inner product is conjugate linear in `z`.
"""->
function inner(x::AbstractVector, y::AbstractVector, z::AbstractVector)
	trapezquad( x, y.*conj(z) )
end

@doc """
	l2norm(x::Vector, y::Vector) -> Number

Approximate the L2 norm squared of `y` over the the interval `[min(x), max(x)]`, using the divisions in `x`.
"""->
function l2norm(x::AbstractVector, y::AbstractVector)
	trapezquad( x, abs2(y) )
end

