"""
    eigval1(A::Matrix) -> Vector

If `A` has 1 as an eigenvalue of multiplicity 1, the corresponding eigenvector is returned.
Otherwise, an error is thrown.
"""
function eigval1(A::AbstractMatrix)::Vector{Float64}
    eigen_values = LinearAlgebra.eigvals(A)

    if !isreal(eigen_values)
        throw(DomainError(eigen_values, "Eigen values should be real"))
    end

    eigen_value_one = findall(eigen_values .≈ 1)

    if length(eigen_value_one) != 1
        throw(DomainError(eigen_value_one, "The eigenvalue 1 should have multiplicity 1"))
    end

    return LinearAlgebra.eigvecs(A)[:, eigen_value_one[1]]
end

