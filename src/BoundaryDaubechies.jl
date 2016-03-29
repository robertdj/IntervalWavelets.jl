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
		coef_mat[row,:] = sqrt2*coef[1:vm]
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
	# TODO: How should E be scaled?

	return E
end


