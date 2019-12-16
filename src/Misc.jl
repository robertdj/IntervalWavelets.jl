"""
    eigval1(A::Matrix) -> Vector

If `A` has 1 as an eigenvalue of multiplicity 1, the corresponding eigenvector is returned.
Otherwise, an error is thrown.
"""
function eigval1(A::AbstractMatrix)
    eigen_values = LinearAlgebra.eigvals(A)

    # TODO: Does this ensure that output type is a real matrix?
    if !isreal(eigen_values)
        throw(DomainError(eigen_values, "Eigen values should be real"))
    end

    eigen_value_one = findall(eigen_values .≈ 1)

    # TODO: Are these good exceptions?
    if isempty(eigen_value_one)
        throw(InvalidStateException)
    end

    if length(eigen_value_one) > 1
        throw(InvalidStateException)
    end

    return LinearAlgebra.eigvecs(A)[:, eigen_value_one[1]]
end

