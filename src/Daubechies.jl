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
	# The min_index bound of the support
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
		return DyadicRationalsVector(0, [1.0; 0.0])
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
	#= @show L = length(y) + (max_index(y) - min_index(y))*2^y.resolution =#
	y_indexvalues = indexvalues(y)
	y_min_index = Int(y_indexvalues[1])
	y_max_index = Int(y_indexvalues[end])

	yy_length = 1 + (y_max_index - y_min_index) * 2^(resolution(y)+1)

	yy = OffsetVector(zeros(yy_length), 2*min_index(y):2*max_index(y))

	yindex = min_index(y):max_index(y)

	i = 2*min_index(y) - 1
	while i < 2*max_index(y)
		# Even indices
		i += 1
		yy[i] = y[div(i, 2)]
		
		# Odd indices
		i += 1
		for k in support(h)
			j = i - k * 2^resolution(y)
			if checkindex(Bool, yindex, j)
				yy[i] += sqrt2 * h[k] * y[j]
			end
		end
	end

	DyadicRationalsVector(y.resolution+1, yy)
end

"""
Compute the scaling function defined by the filter `h` at the dyadic
rationals of resolution `R`.
"""
function DaubScaling(h::InteriorFilter, R::Int)
	phi = DaubScaling(h)
	for res in 1:R
		phi = DaubScaling(phi, h)
	end

	return phi
end

