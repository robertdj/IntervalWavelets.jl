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
		vm = WT.vanishingmoments( WT.eval(parse(lowername)) )
		C = scalingfilters(vm)
	end
end

