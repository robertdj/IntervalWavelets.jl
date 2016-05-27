@doc """
	boundary_coef_mat(F::BoundaryFilter) -> Matrix

The boundary coefficients collected in a matrix where the `i`'th row
contains the coefficients of the `i`'th boundary scaling function.
"""->
function boundary_coef_mat(F::BoundaryFilter)
	const vm = van_moment(F)
	coef_mat = zeros(Float64, vm, vm)

	for row = 1:vm
		coef = bfilter(F, row-1)
		@inbounds coef_mat[row,:] = sqrt2*coef[1:vm]
	end

	return coef_mat
end

@doc """
	DaubScaling(p::Integer, side::Char, R::Integer) -> x, Y

Compute the boundary scaling function at resolution `R` on `side` with `p` vanishing moments.
"""->
function DaubScaling(p::Integer, side::Char, R::Integer)
	B = bfilter(p, side)
	IF = ifilter(p, true)

	x = dyadic_rationals( support(B), R )
	Y = DaubScaling(B, IF, R)

	return x, Y'
end

@doc """
	DaubScaling(B::BoundaryFilter) -> Vector

Compute the boundary scaling function values at 0.
"""->
function DaubScaling(B::BoundaryFilter)
	boundary_mat = boundary_coef_mat(B)
	E = eigval1(boundary_mat)
	# TODO: How to scale E?

	return E
end

@doc """
	DaubScaling(B, I) -> Matrix

Compute the boundary scaling function defined by boundary filter `B` and
interior filter `F` values at the non-zero integers in their support.

The ouput is a matrix where the `k`'th row are the functions values of the `k-1` scaling function.
"""->
function DaubScaling(B::BoundaryFilter, IF::InteriorFilter)
	const internal = DaubScaling(IF)
	const IS = support(IF)
	const BS = support(B)

	# Function values at 0
	const vm = van_moment(B)
	const xvals = integers(B)
	Y = zeros(Float64, vm, length(xvals)+1)
	Y[:,x2index(0,BS)] = DaubScaling(B)

	const xvals = integers(B)
	# The translations of the interior scaling function differ for the two sides
	const interior_start = ( side(B) == 'L' ? vm : -vm-1 )
	const interior_iter = ( side(B) == 'L' ? 1 : -1 )

	for x in xvals
		xindex = x2index(x, BS)
		doublex = 2*x
		doublex_index = x2index(doublex, BS)

		for k in 1:vm
			# TODO: copy! ?
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
					Y[k,xindex] += sqrt2 * filterk[fidx] * internal[ x2index(interior_arg,IS) ]
				end
				interior_arg -= interior_iter
			end
		end
	end

	return Y
end

@doc """
	DaubScaling(B, I, R) -> Matrix

Compute the boundary scaling function defined by boundary filter `B` and interiorfilter `I` values at the dyadic rationals up to resolution `R`in their support.

The ouput is a matrix where the `k`'th row are the functions values of the `k-1` scaling function.
"""->
function DaubScaling(B::BoundaryFilter, IF::InteriorFilter, R::Int)
	@assert R >= 0

	const internal = DaubScaling(IF,R)
	const Ny = length(internal)
	const vm = van_moment(B)
	Y = zeros(Float64, vm, Ny)

	const IS = support(IF)
	const BS = support(B)

	# Base level
	cur_idx = dyadic_rationals(BS, R, 0)
	Y[:,cur_idx] = DaubScaling(B, IF)

	const interior_start = ( side(B) == 'L' ? vm : -vm-1 )
	const interior_iter = ( side(B) == 'L' ? 1 : -1 )

	# Recursion: Fill remaining levels
	const xvals = dyadic_rationals(BS, R)
	for L in 1:R
		# Indices of x values on scale L
		cur_idx = dyadic_rationals(BS, R, L)

		for xindex in cur_idx
			doublex = 2*xvals[xindex]
			doublex_index = x2index(doublex, BS, R)

			for k in 1:vm
				# TODO: copy! ?
				filterk = bfilter(B, k-1)

				# Boundary contribution
				# TODO: The support depends on k
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

