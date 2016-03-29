@doc """
The support of a Daubechies scaling function is an interval with integer end points.
"""->
type DaubSupport
	left::Integer
	right::Integer

	DaubSupport(left,right) = right <= left ? error("Not an interval") : new(left, right)
end

left(I::DaubSupport) = I.left
right(I::DaubSupport) = I.right
Base.length(I::DaubSupport) = right(I) - left(I)

@doc """
	support(C) -> DaubSupport

Return the length of the support of the scaling function defined by the
filter vector `C`.
"""->
function support(C::Vector{Float64})
	return DaubSupport(0, length(C)-1)
end

@doc """
	isinside(x, I::DaubSupport) -> Bool

Returns `true` if `x` is inside `I`.
"""->
isinside(x, I::DaubSupport) = left(I) <= x <= right(I)

# Convert between values and indices of a vector with the integers in
# the support I
x2index(x::Integer, I::DaubSupport) = x + 1 - left(I)
index2x(idx::Integer, I::DaubSupport) = idx - 1 + left(I)

# Convert between values and indices of a vector with dyadic rationals
# at resolution R in the support I
x2index(x, I::DaubSupport, R::Integer) = Int( (x-left(I))*2^R ) + 1
index2x(idx::Integer, I::DaubSupport, R::Integer) = (idx - 1)/2^R + left(I)

