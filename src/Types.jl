typealias DaubSupport Tuple{Int,Int}
left(I::DaubSupport) = I[1]
right(I::DaubSupport) = I[2]

#= Base.length(I) = right(I) - left(I) =#
@doc """
	support(C) -> DaubSupport

Return the length of the support of the scaling function defined by the
filter vector `C`.
"""->
function support(C::Vector{Float64})
	return (0,length(C) - 1)
end


isinside(x, I::DaubSupport) = left(I) <= x <= right(I)


# Convert between values and indices of a vector with the integers in
# the support I
x2index(x::Integer, I::DaubSupport) = x + 1 - left(I)
index2x(idx::Integer, I::DaubSupport) = idx - 1 + left(I)

# Convert between values and indices of a vector with dyadic rationals
# at resolution R in the support I
x2index(x, I::DaubSupport, R::Integer) = Int( (x-left(I))*2^R ) + 1
index2x(idx::Integer, I::DaubSupport, R::Integer) = (idx - 1)/2^R + left(I)

