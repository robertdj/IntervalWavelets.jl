@doc """
	eigval1(A::Matrix) -> Vector

Returns a unit eigenvector of `A` for the eigenvalue 1, if this eigenspace is 1D.
Otherwise, return an error.
"""->
function eigval1(A::Matrix)
	evals = eigvals(A)

	# Find index of eigenvalue 1
	isevals1 = map( x->isapprox(x,1.0), evals )
	eval1_index = find(isevals1)

	if isempty(eval1_index)
		error("1 is not an eigenvalue")
	elseif length(eval1_index) >= 2
		error("The eigenspace of 1 is more than 1D")
	end

	# Eigenvector of 1
	D = evals[eval1_index[]]*eye(A)
	evec = nullspace( A - D )

	return evec
end

