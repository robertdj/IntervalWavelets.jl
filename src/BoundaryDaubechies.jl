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
	eigval1(A::Matrix) -> Vector

Returns a unit eigenvector of `A` for the eigenvalue 1.
It is assumed that the eigenvector is unique.
"""->
function eigval1(A::Matrix)
	eval, evec = eig(A)

	# Find index of eigenvalue 1
	idx = 1
	for v in eval
		if isapprox(v, 1.0)
			break
		end
		idx += 1
	end

	if idx > length(eval)
		error("1 is not an eigenvalue")
	end

	return evec[:,idx]
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


