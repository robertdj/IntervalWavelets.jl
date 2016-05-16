@doc """
	eigval1(A::Matrix) -> Vector

Returns a unit eigenvector of `A` for the eigenvalue 1 if this eigenspace is 1D.
Otherwise, return an error.
"""->
function eigval1(A::Matrix)
	evals, evecs = eig(A)

	# Find index of eigenvalue 1
	isevals1 = map( x->isapprox(x,1.0), evals )
	eval1_index = find(isevals1)

	if isempty(eval1_index)
		error("1 is not an eigenvalue")
	elseif length(eval1_index) >= 2
		error("The eigenspace of 1 is more than 1D")
	end

	return evecs[:,eval1_index[]]
end


@doc """
	waveparse(wavename::AbstractString) -> Vector

Return `wavename`'s filter.

	waveparse(wavename::AbstractString, true) -> ScalingFilters

Return `wavename`'s internal and boundary filters.

Currently, only Daubechies wavelets are supported and they are denoted 

- `"Haar"`.
- `"dbN"` where `N` is an integer between 1 and 8.
"""->
function waveparse(wavename::AbstractString, boundary::Bool=false)
	# TODO: Enforce Daubechies only
	lowername = lowercase( wavename )

	# The Haar wavelet does not have boundary issues
	if lowername == "haar" || lowername == "db1"
		boundary = false
	end

	if !boundary
		C = wavelet( WT.eval(parse(lowername)) ).qmf
	else
		vm = van_moment(lowername)
		C = scalingfilters(vm)
	end
end

function van_moment(wavename::AbstractString)
	lowername = lowercase( wavename )
	WT.vanishingmoments( WT.eval(parse(lowername)) )
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
x2index(x, S::DaubSupport, R::Integer) = Int( (x-left(S))*2^R ) + 1
index2x(idx::Integer, S::DaubSupport, R::Integer) = (idx - 1)*2.0^(-R) + left(S)

