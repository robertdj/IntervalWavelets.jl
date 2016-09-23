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
	DaubScaling(p, R) -> x, y

A Daubechies `p` scaling function evaluated in the dyadic rationals at resolution `R`.
"""->
function DaubScaling(p::Integer, R::Integer, symmlet::Bool=true)
	IF = ifilter(p, symmlet)
	supp = support(IF)
	x = dyadic_rationals(supp, R)
	phi = DaubScaling( IF, R )

	return x, phi
end

@doc """
	dyadic_dil_matrix(C::Vector) -> Matrix

The "dyadic dilation matrix" `D` of the filter `C`: 
`D[i,j] = C[2i-j]`.
"""->
function dyadic_dil_matrix(C::Vector{Float64})
	NC = length(C)
	sz = NC - 2
	dydil_mat = zeros(Float64, sz, sz)

	for nj in 1:sz, ni in 1:sz
		Cidx = 2*ni - nj + 1
		if 1 <= Cidx <= NC
			dydil_mat[ni, nj] = sqrt2*C[Cidx]
		end
	end

	return dydil_mat
end


#=
	DaubScaling(IF::InteriorFilter) -> Vector

Compute function values of the scaling function defined by the filter `IF` at the integers in the support.
=#
function DaubScaling(IF::InteriorFilter)
	if van_moment(IF) == 1
		return [1.0; 0.0]
	end

	C = coef(IF)
	L = dyadic_dil_matrix( C )

	# Eigenvector of eigenvalue 1: 
	E = zeros(Float64, length(C))
	E[2:end-1] = eigval1(L)

	# Normalize scaling function in L2
	scale!(E, 1/sum(E))

	return E
end

@doc """
	DaubScaling(C::InteriorFilter, R::Int) -> Vector

Compute function values of the scaling function defined by the filter
`C` at the dyadic rationals of resolution `R` in the support.
"""->
function DaubScaling(IF::InteriorFilter, R::Integer)
	supp = support(IF)
	# There are 2^R points on each unit + endpoint
	Nx = length(supp)*2^R + 1
	phi = zeros(Float64, Nx)

	# Base level
	cur_idx = dyadic_rationals(supp, R, 0)
	phi[cur_idx] = DaubScaling(IF)

	# Recursion: Fill remaining levels
	coeff = coef(IF)
	Lsupp = left(supp)
	Rsupp = right(supp)
	for L in 1:R
		# Indices of x values on scale L
		cur_idx = dyadic_rationals(supp, R, L)

		for xindex in cur_idx
			twox = 2*index2x( xindex, supp, R )

			for k in Lsupp:Rsupp
				if isinside(twox-k, supp)
					twox_index = x2index( twox-k, supp, R )
					phi[xindex] += sqrt2 * coeff[k-Lsupp+1] * phi[twox_index]
				end
			end
		end
	end

	return phi
end

