"""
	boundary_coef_mat(F::BoundaryFilter) -> Matrix

The boundary coefficients collected in a matrix where the `i`'th row
contains the coefficients of the `i`'th boundary scaling function.
"""
function boundary_coef_mat(F::BoundaryFilter)
	vm = van_moment(F)
	coef_mat = zeros(Float64, vm, vm)

	for row in 1:vm
		coef = filter(F, row-1)
		for col in 1:vm
			coef_mat[row,col] = sqrt2*coef[col]
		end
	end

	return coef_mat
end

"""
	DaubScaling(p::Integer, side::Char, R::Integer) -> x, Y

Compute the boundary scaling function at resolution `R` on `side` with `p` vanishing moments.
"""
function DaubScaling(p::Integer, side::Char, R::Integer)
	B = bfilter(p, side)
	IF = ifilter(p, true)

	x = dyadic_rationals( support(B), R )
	Y = DaubScaling(B, IF, R)

	return x, Y'
end

#=
	DaubScaling(B::BoundaryFilter) -> Vector

Compute the boundary scaling function values at 0.
=#
function DaubScaling(B::BoundaryFilter)
	boundary_mat = boundary_coef_mat(B)
	E = eigval1(boundary_mat)
	# TODO: How to scale E?

	return E
end

function DaubScaling(H::BoundaryFilter, h::InteriorFilter)
	if (p = van_moment(h)) != van_moment(H)
		throw(AssertionError())
	end

	# The DyadicRationalsVector is not an AbstractArray and cannot be
	# used in an OffsetVector
	Y = Vector{DyadicRationalsVector}(p)
	Y_indices = OffsetVector{UnitRange{Int64}}(0:p-1)
	for k in 0:p-1
		Y_indices[k] = 0:p+k
		Y[k+1] = DyadicRationalsVector(0, OffsetVector(zeros(p+k+1), 0:p+k))

		# TODO: How to compute function values at 0?
		Y[k+1][0] = NaN
	end

	# Interior scaling function
	y = DaubScaling(h)
	y_indices = linearindices(y)

	# Loop over the integers with non-zero function values, which
	# depends on the function index (k)
	for x in 2*p-2:-1:1
		for k in p-1:-1:max(x-p+1,0)
			yval = 0.0

			# Boundary contribution
			for Hi in 0:p-1
				Yi = 2*x
				if checkindex(Bool, Y_indices[Hi], Yi)
					yval += sqrt2 * filter(H, k)[Hi] * Y[Hi+1][Yi]
				end
			end

			# Interior contribution
			for Hi in p:p+2*k
				yi = 2*x - Hi
				if checkindex(Bool, y_indices, yi)
					yval += sqrt2 * filter(H, k)[Hi] * y[yi]
				end
			end

			Y[k+1][x] = yval
		end
	end

	return Y
end


function DaubScaling(H::BoundaryFilter, h::InteriorFilter, R::Int)
	if R < 0
		throw(DomainError())
	end

	Y = DaubScaling(H, h)
	y = DaubScaling(h)
	for res in 1:R
		Y = DaubScaling(Y, y, H, h)
		y = DaubScaling(y, h)
	end

	return Y
end


#=
	DaubScaling(B, I) -> Matrix

Compute the boundary scaling function defined by boundary filter `B` and
interior filter `I` at the non-zero integers in their support.

The ouput is a matrix where the `k`'th row are the functions values of the `k-1` scaling function.
=#
function DaubScaling1(B::BoundaryFilter, IF::InteriorFilter)
	internal = DaubScaling(IF)
	IS = support(IF)
	BS = support(B)

	# Function values at 0
	vm = van_moment(B)
	xvals = integers(B)
	Y = zeros(Float64, vm, length(xvals)+1)
	Y[:, x2index(0,BS)] = DaubScaling(B)

	xvals = integers(B)
	# The translations of the interior scaling function differ for the two sides
	interior_start = ( side(B) == 'L' ? vm : -vm-1 )
	interior_iter = ( side(B) == 'L' ? 1 : -1 )

	for x in xvals
		xindex = x2index(x, BS)
		doublex = 2*x
		doublex_index = x2index(doublex, BS)

		for k in 1:vm
			filterk = bfilter(B, k-1)

			# Boundary contribution
			if isinside(doublex, BS)
				for l in 1:vm
					Y[k, xindex] += sqrt2 * filterk[l] * Y[l, doublex_index]
				end
			end

			# Interior contribution
			fidx = vm
			interior_arg = doublex - interior_start
			flength = vm + 2*(k-1) #length(filterk)
			while fidx <= flength
				fidx += 1
				if isinside(interior_arg, IS)
					Y[k,xindex] += sqrt2 * filterk[fidx] * internal[ x2index(interior_arg,IS) ]
				end
				interior_arg -= interior_iter
			end
		end
	end

	return Y
end

"""
	DaubScaling(B, I, R) -> Matrix

Compute the boundary scaling function defined by boundary filter `B` and interior filter `I` at the dyadic rationals up to resolution `R` in their support.

The ouput is a matrix where the `k`'th row are the functions values of the `k-1` scaling function.
"""
function DaubScaling1(B::BoundaryFilter, IF::InteriorFilter, R::Int)
	R >= 0 || throw(DomainError())

	internal = DaubScaling(IF,R)
	Ny = length(internal)
	vm = van_moment(B)
	Y = zeros(Float64, vm, Ny)

	IS = support(IF)
	BS = support(B)

	# Base level
	cur_idx = dyadic_rationals(BS, R, 0)
	Y[:,cur_idx] = DaubScaling(B, IF)

	interior_start = ( side(B) == 'L' ? vm : -vm-1 )
	interior_iter = ( side(B) == 'L' ? 1 : -1 )

	# Recursion: Fill remaining levels
	xvals = dyadic_rationals(BS, R)
	for L in 1:R
		# Indices of x values on scale L
		cur_idx = dyadic_rationals(BS, R, L)

		for xindex in cur_idx
			doublex = 2*xvals[xindex]
			doublex_index = x2index(doublex, BS, R)

			for k in 1:vm
				filterk = bfilter(B, k-1)

				# Boundary contribution
				if isinside(doublex, BS)
					for l in 1:vm
						Y[k,xindex] += sqrt2 * filterk[l] * Y[l,doublex_index]
					end
				end

				# Interior contribution
				fidx = vm
				interior_arg = doublex - interior_start
				flength = vm + 2*(k-1) #length(filterk)
				while fidx <= flength
					fidx += 1
					if isinside(interior_arg, IS)
						Y[k,xindex] += sqrt2 * filterk[fidx] * internal[ x2index(interior_arg,IS,R) ]
					end
					interior_arg -= interior_iter
				end
			end
		end
	end

	# Hack until the right 0 values are computed
	# TODO: Fix this
	if side(B) == 'L'
		Y[:,1] = Y[:,2]
	else
		Y[:,end] = Y[:,end-1]
	end

	return Y
end

