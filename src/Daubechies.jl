# ------------------------------------------------------------
# Scaling functions

"""
	HaarScaling( xi[, J, k] ) -> Float

The Haar scaling function evaluated in `xi` at level `J` and translation `k`.
By default, `J=0` and `k=0`.
"""
function HaarScaling(xi)
	zero(xi) <= xi < one(xi) ? 1.0 : 0.0
end

"""
	DaubScaling(p, R) -> x, y

A Daubechies `p` scaling function evaluated in the dyadic rationals at resolution `R`.
"""
function DaubScaling(p::Integer, R::Integer, symmlet::Bool=true)
	IF = ifilter(p, symmlet)
	supp = support(IF)
	x = dyadic_rationals(supp, R)
	phi = DaubScaling( IF, R )

	return x, phi
end

"""
	dyadic_dilation_matrix(h::InteriorFilter)

Compute the matrix whose eigenvalues is the wavelet values at the
integers in its support.
See the `doc` folder.
"""
function dyadic_dilation_matrix(h::InteriorFilter)
	sz = length(h) - 2
	ddmat = zeros(Float64, sz, sz)
	# The lower bound of the support
	l = support(h)[1]

	for j in 1:sz, i in 1:sz
		ddmat[i, j] = sqrt2*h[2*i - j + l]
	end

	return ddmat
end


#=
	DaubScaling(IF::InteriorFilter) -> Vector

Compute function values of the scaling function defined by the filter `IF` at the integers in the support.
=#
function DaubScaling(h::InteriorFilter)
	# Haar wavelet. Not covered in dyadic_dilation_matrix because there
	# are no integers in the interior of the support
	if van_moment(h) == 1
		return [1.0; 0.0]
	end

	# Non-zero function values
	support_idx = support(h)
	E = OffsetVector(zeros(Float64, length(h)), support_idx)

	nonzero_idx = support_idx[2]:support_idx[end-1]
	E[nonzero_idx] = h |> 
		dyadic_dilation_matrix |> 
		eigval1 |> 
		x -> scale!(x, 1/sum(x)) # Normalize scaling function in L2

	return DyadicRationalsVector(0, E)
end

"""
	DaubScaling(C::InteriorFilter, R::Int) -> Vector

Compute function values of the scaling function defined by the filter
`C` at the dyadic rationals of resolution `R` in the support.
"""
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

