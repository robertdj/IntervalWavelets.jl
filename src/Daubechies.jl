# ------------------------------------------------------------
# Scaling functions

@doc """
	HaarScaling( xi[, J, k] ) -> Float

The Haar scaling function evaluated in `xi` at level `J` and translation `k`.
By default, `J=0` and `k=0`.
"""->
function HaarScaling(xi::Real)
	zero(xi) <= xi < one(xi) ? 1.0 : 0.0
end

@doc """
	DaubScaling(N, R) -> Vector, Vector

A Daubechies `N` scaling function evaluated in the dyadic rationals at resolution `R`.
"""->
function DaubScaling(N::Int, R::Int)
	const C = wavelet( WT.Daubechies{N}() ).qmf
	const supp = support(C)
	x = dyadic_rationals(supp, R)
	phi = DaubScaling( C, R )

	return x, phi
end

@doc """
	dyadic_dil_matrix(C::Vector) -> Matrix

The "dyadic dilation matrix" `D` of the filter `C` with `D[i,j] = C[2i-j]`.
"""->
function dyadic_dil_matrix(C::Vector{Float64})
	const NC = length(C)
	dydil_mat = zeros(Float64, NC, NC)

	const sqrt2 = sqrt(2)
	for nj = 1:NC, ni = 1:NC
		Cidx = 2*ni - nj
		# TODO: Avoid this check?
		if 1 <= Cidx <= NC
			dydil_mat[ni, nj] = sqrt2*C[Cidx]
		end
	end

	return dydil_mat
end


@doc """
	DaubScaling(C::Vector) -> Vector

Compute function values of the scaling function defined by the filter `C` at the integers in the support.
"""->
function DaubScaling(C::Vector{Float64})
	L = dyadic_dil_matrix(C)

	# Eigenvector of eigenvalue 1
	E = nullspace( L-eye(length(C)) )

	# The scaling function should be normalized in L2
	scale!(E, 1/sum(E))

	# The first and last entry are both 0
	# TODO: Not for Haar
	# TODO: Don't compute them
	E[1] = E[end] = 0.0

	return E
end

@doc """
	DaubScaling(C::Vector, R::Int) -> Vector

Compute function values of the scaling function defined by the filter
`C` at the dyadic rationals of resolution `R` in the support.
"""->
function DaubScaling(C::Vector{Float64}, R::Int)
	const supp = support(C)
	const x = dyadic_rationals(supp, R)
	const Nx = length(x)
	phi = zeros(Float64, Nx)

	# Base level
	cur_idx = dyadic_rationals( supp, R, 0)
	phi[cur_idx] = DaubScaling(C)

	# Recursion: Fill remaining levels
	const NC = length(C)
	const sqrt2 = sqrt(2)
	for L = 1:R
		# Indices of x values on scale L
		cur_idx = dyadic_rationals(supp, R, L)

		for phin = cur_idx
			for Ck = 1:NC
				pidx = dyadic_parent(phin, Ck-1, R)
				# TODO: Calculate only the necessary pidx
				if 1 <= pidx <= Nx
					phi[phin] += sqrt2*C[Ck]*phi[pidx]
				end
			end
		end
	end

	return phi
end

