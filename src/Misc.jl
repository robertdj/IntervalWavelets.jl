function unit(x::AbstractVector, entry::Integer)
	unit(eltype(x), length(x), entry)
end

function unit(length::Integer, entry::Integer)
	unit(Float64, length, entry)
end

function unit(T::Type, length::Integer, entry::Integer)
	u = zeros(T, length)
	u[entry] = one(T)
	return u
end

"""
	eigval1(A::Matrix) -> Vector

Returns a unit eigenvector of `A` for the eigenvalue 1 if this eigenspace is 1D.
Otherwise, return an error.
"""
function eigval1(A::DenseMatrix{Float64})
	E = eigfact(A)
	# Types are hardcoded; otherwise instability occurs
	isreal(E[:values]) || throw(DomainError())
	vals = E[:values]::Vector{Float64}
	vecs = E[:vectors]::Matrix{Float64}

	# Find index of eigenvalue 1
	eval1_index = Vector{Int64}()
	for i in 1:length(vals)
		if isapprox(vals[i], 1.0)
			push!(eval1_index, i)
		end
	end

	if isempty(eval1_index)
		error("1 is not an eigenvalue")
	elseif length(eval1_index) >= 2
		error("The eigenspace of 1 is more than 1D")
	end

	return vecs[:,eval1_index[]]
end

function van_moment(wavename::String)
	lowername = lowercase( wavename )
	WT.vanishingmoments( WT.eval(parse(lowername)) )::Int64
end


"""
	support(wavename) -> DaubSupport

The support of the Daubechies scaling function `wavename`.
"""
function support(wavename::String)
	vm = van_moment(wavename)
	return DaubSupport(-vm+1, vm)
end


"""
	isinside(x, S::DaubSupport) -> Bool

Returns `true` if `x` is inside `S`.
"""
@inline isinside(x, S::DaubSupport) = left(S) <= x <= right(S)

"""
Convert between values and indices of a vector with the integers in the support S
"""
@inline x2index(x::Integer, S::DaubSupport) = x + 1 - left(S)
@inline index2x(idx::Integer, S::DaubSupport) = idx - 1 + left(S)

"""
Convert between values and indices of a vector with dyadic rationals
at resolution R in the support S
"""
@inline x2index(x, S::DaubSupport, R::Integer) = Int( (x-left(S))*2^R ) + 1
@inline index2x(idx::Integer, S::DaubSupport, R::Integer) = (idx - 1)*2.0^(-R) + left(S)

