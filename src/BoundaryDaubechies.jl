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
			coef_mat[row, col] = sqrt2*coef[col]
		end
	end

	return coef_mat
end

"""
	DaubScaling(p::Integer, side::Char, R::Integer) -> x, Y

Compute the boundary scaling function at resolution `R` on `side` with `p` vanishing moments.
"""
function DaubScaling(p::Integer, side::Char, R::Integer, phase::String="symmlet")
	B = bfilter(p, side)
	IF = ifilter(p, phase)

	x = dyadic_rationals( support(B), R )
	Y = DaubScaling(B, IF, R)

	return x, Y'
end

# TODO: Don't keep this function if we don't use it
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
	# used in an OffsetVector. In the remainder we have to +1 when indexing Y
	Y = Vector{DyadicRationalsVector}(p)
	Y_indices = OffsetVector{UnitRange{Int64}}(0:p-1)
	for k in 0:p-1
		Y_indices[k] = support_integers(H, k)
		Y[k+1] = DyadicRationalsVector(0, OffsetVector(zeros(p+k+1), 
													   Y_indices[k]))

		# TODO: How to compute function values at 0?
		Y[k+1][0] = NaN
	end

	# Interior scaling function
	y = DaubScaling(h)
	y_indices = linearindices(y)

	# Loop over the integers with non-zero function values, which
	# depends on the function index (k)
	Y_support = loop_support_integers(H)

	for x in Y_support
		for k in p-1:-1:0
			if !checkindex(Bool, Y_indices[k], x)
				continue
			end

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
				yi = interior_index(H, x, Hi)
				if checkindex(Bool, y_indices, yi)
					yval += sqrt2 * filter(H, k)[Hi] * y[yi]
				end
			end

			Y[k+1][x] = yval
		end
	end

	return Y
end


function DaubScaling(Y::Vector{DyadicRationalsVector}, y::DyadicRationalsVector,
				     H::BoundaryFilter, h::InteriorFilter)
	if (p = van_moment(h)) != van_moment(H)
		throw(AssertionError())
	end
	# TODO: Check all resolutions in Y
	if (R = resolution(y)) != resolution(Y[1])
		throw(AssertionError())
	end

	y_indices = linearindices(y)
	#= Y_indices = map(linearindices, Y) =#
	#= Y_indices = OffsetVector(map(linearindices, Y), 0:p-1) =#
	Y_indices = OffsetVector{UnitRange{Int64}}(0:p-1)
	for k in 0:p-1
		Y_indices[k] = linearindices(Y[k+1])
	end

	unitstep = 2^R
	Z = Vector{DyadicRationalsVector}(p)

	for k in 0:p-1
		Y_support = support(Y[k+1])
		Y_min_support = Int(Y_support[1])
		Y_max_support = Int(Y_support[end])
		#= y_indices = linearindices(Y[k+1]) =#

		z_length = 1 + (Y_max_support - Y_min_support) * 2^(R+1)
		z_min_index = 2*Y_indices[k][1]
		z_max_index = 2*Y_indices[k][end]
		z = OffsetVector(Array{Float64}(z_length), z_min_index:z_max_index)

		# The even indices are inherited from the input. The odd indices are
		# computed using the dilation equation
		# TODO: In Julia v0.7, use the `first` index
		zi = z_min_index
		z[zi] = Y[k+1][Y_indices[k][1]]
		while zi < z_max_index
			# ------------------------------------------------------------------
			# Odd indices
			zi += 1
			zval = 0.0

			# Boundary contribution
			for Hi in 0:p-1
				Yi = zi
				if checkindex(Bool, Y_indices[Hi], Yi)
					zval += sqrt2 * filter(H, k)[Hi] * Y[Hi+1][Yi]
				end
			end

			# Interior contribution
			for Hi in p:p+2*k
				yi = interior_index(H, zi, Hi, unitstep)
				if checkindex(Bool, y_indices, yi)
					zval += sqrt2 * filter(H, k)[Hi] * y[yi]
				end
			end

			z[zi] = zval

			# ------------------------------------------------------------------
			# Even indices
			zi += 1
			z[zi] = Y[k+1][div(zi, 2)]
		end

		Z[k+1] = DyadicRationalsVector(R+1, z)
	end

	return Z
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



# ------------------------------------------------------------------------------

function support_integers(H::BoundaryFilter, k::Int)
	p = van_moment(H)
	if side(H) == 'L'
		return 0:p+k
	else
	#= elseif side(H) == 'R' =#
		# Type unstable with elseif as the compiler thinks that "Void"
		# is a possible return type. Change if this function is
		# @generated
		return -p-k:0
	end
end

function loop_support_integers(H::BoundaryFilter)
	p = van_moment(H)
	if side(H) == 'L'
		return 2*p-2:-1:1
	else
		return -2*p+2:-1
	end
end

function interior_index(H::BoundaryFilter, x::Int, Hi::Int)
	if side(H) == 'L'
		return 2*x - Hi
	else
		return 2*x + Hi + 1
	end
end

function interior_index(H::BoundaryFilter, zi::Int, Hi::Int, unitstep::Int)
	if side(H) == 'L'
		return zi - Hi * unitstep
	else
		return zi + (Hi + 1) * unitstep
	end
end

