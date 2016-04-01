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
	DaubScaling(B::BoundaryFilter) -> Vector

Compute the boundary scaling function values at 0.
"""->
function DaubScaling(B::BoundaryFilter)
	boundary_mat = boundary_coef_mat(B)
	E = eigval1(boundary_mat)
	# TODO: How to scale E?

	return E
end


function DS(x, k::Int, F::ScalingFilters)
	if x == 0
		return 0
	end

	R = dyadic_level(x)

	vm = van_moment(F)
	@assert 0 <= k < vm
	IS, BS = support(F)

	phi = DaubScaling(F.internal, R)
	xx = dyadic_rationals(IS, R)

	if isinside(x, BS)
		Y = zeros(Float64, length(bfilter(F.left,k)))
		for l = 0:vm-1
			Y[l+1] = DS(2*x, l, F)
		end
		for m = vm:(vm+2*k)
			if isinside(2*x-m, IS)
				Y[m+1] = phi[ x2index(2*x-m,IS,R) ]
			end
		end

		@show x, k, Y
		println( bfilter(F.left,k) )
		#= @bp =#
		return sqrt(2)*dot(Y, bfilter(F.left,k))
	else
		return zero(float(x))
	end
end

function DS(x::Vector, k::Int, F::ScalingFilters)
	map(x -> DS(x,k,F), x)
end

