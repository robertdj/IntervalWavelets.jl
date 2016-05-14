@doc """
The support of a Daubechies scaling function is an interval.
"""->
type DaubSupport
	left::Integer
	right::Integer

	DaubSupport(left,right) = right <= left ? error("Not an interval") : new(left, right)
end

# TODO: WARNING: imported binding for I overwritten in module Main
left(S::DaubSupport) = S.left
right(S::DaubSupport) = S.right
Base.length(S::DaubSupport) = right(S) - left(S)

@doc """
	support(C) -> DaubSupport

The support of the Daubechies scaling function defined by the filter vector `C`.
"""->
function support(C::Vector{Float64})
	#= vm = div(length(C)+1, 2) =#
	vm = length(C) / 2
	return DaubSupport(-vm+1, vm)
end

@doc """
	support(wavename) -> DaubSupport

The support of the Daubechies scaling function `wavename`.
"""->
function support(wavename::AbstractString)
	vm = van_moment(wavename)
	return DaubSupport(-vm+1, vm)
end

@doc """
Compute support of the scaling function/wavelet at scale `J` and with translation `k` from the support `S` of the father/mother.
"""->
function support(S::DaubSupport, J::Integer, k::Integer)
	L = 2.0^(-J)*(left(S) + k)
	R = 2.0^(-J)*(right(S) + k)
	DaubSupport(L, R)
end

@doc """
	isinside(x, S::DaubSupport) -> Bool

Returns `true` if `x` is inside `S`.
"""->
isinside(x, S::DaubSupport) = left(S) <= x <= right(S)

# Convert between values and indices of a vector with the integers in the support S
x2index(x::Integer, S::DaubSupport) = x + 1 - left(S)
index2x(idx::Integer, S::DaubSupport) = idx - 1 + left(S)

# Convert between values and indices of a vector with dyadic rationals
# at resolution R in the support S
x2index(x, S::DaubSupport, R::Integer) = Int( (x-left(S)) )*2^R + 1
index2x(idx::Integer, S::DaubSupport, R::Integer) = (idx - 1)*2.0^(-R) + left(S)

