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
function DaubScaling(p::Integer, R::Integer, phase::String="symmlet")
	h = ifilter(p, phase)
	DaubScaling(h, R)
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


"""
	DaubScaling(h::InteriorFilter)

Compute function values of the scaling function defined by the filter `h` at 
the integers in the support.
"""
function DaubScaling(h::InteriorFilter)
	# Haar wavelet. Not covered in dyadic_dilation_matrix because there
	# are no integers in the interior of the support
	if van_moment(h) == 1
		return DyadicRationalsVector(0, OffsetVector([1.0; 0.0], 0:1))
	end

	support_idx = support(h)
	y = OffsetVector(zeros(Float64, length(h)), support_idx)

	# Non-zero function values
	nonzero_idx = support_idx[2]:support_idx[end-1]
	y[nonzero_idx] = h |> 
		dyadic_dilation_matrix |> 
		eigval1 |> 
		x -> scale!(x, 1/sum(x)) # Normalize scaling function in L2

	return DyadicRationalsVector(0, y)
end

"""
Increase the resolution of the DaubScaling scaling function by one.
"""
function DaubScaling(y::DyadicRationalsVector, h::InteriorFilter)
	y_support = support(y)
	y_min_support = Int(y_support[1])
	y_max_support = Int(y_support[end])
	y_indices = linearindices(y)

	z_length = 1 + (y_max_support - y_min_support) * 2^(resolution(y)+1)
	z_min_index = 2*y_indices[1]
	z_max_index = 2*y_indices[end]
	z = OffsetVector(zeros(z_length), z_min_index:z_max_index)

	# The even indices are inherited from the input. The odd indices are
	# computed using the dilation equation
	# The first value is always zero (due to the continuity and compact
	# support), so we don't copy this (makes the code easier)
	unitstep = 2^resolution(y)
	zi = z_min_index
	while zi < z_max_index
		# Odd indices
		zi += 1
		zval = 0.0
		for hi in support(h)
			yi = zi - hi * unitstep
			if checkindex(Bool, y_indices, yi)
				zval += sqrt2 * h[hi] * y[yi]
			end
		end
		z[zi] = zval

		# Even indices
		zi += 1
		z[zi] = y[div(zi, 2)]
	end

	DyadicRationalsVector(resolution(y)+1, z)
end

"""
Compute the scaling function defined by the filter `h` at the dyadic
rationals of resolution `R`.
"""
function DaubScaling(h::InteriorFilter, R::Int)
	if R < 0
		throw(DomainError())
	end

	phi = DaubScaling(h)
	for res in 1:R
		phi = DaubScaling(phi, h)
	end

	return phi
end

