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

function DS(x::Real, C::Vector{Float64})
	const S = support(C)
	if !isinside(x, S)
		return 0.0
	end

	const N = length(C)
	const p = div(N,2)

	if isinteger(x)
		y = DaubScaling(C)
		return y[Int(x)+p]
	else
		Y = Array{Float64}(N)
		idx = 0
		for k in -p+1:p
			idx += 1
			Y[idx] = DS(2*x-k, C)
		end
		return sqrt2*dot(C,Y)
	end
end

function DS(x::Real, k::Int, F::ScalingFilters)
	const IS, BS = support(F)
	if !isinside(x, BS) || x == 0
		return 0.0
	end

	const vm = van_moment(F)
	@assert 0 <= k < vm

	const L = bfilter(F.left,k)
	Y = zeros(Float64, length(L))

	for l = 0:vm-1
		Y[l+1] = DS(2*x, l, F)
	end

	for m = vm:(vm+2*k)
		Y[m+1] = DS(2*x-m, F.internal)
	end

	return sqrt2*dot(Y, L)
end

function DS(x::Vector, C::Vector{Float64})
	map(x -> DS(x,C), x)
end

function DS(x::Vector, k::Int, F::ScalingFilters)
	map(x -> DS(x,k,F), x)
end

