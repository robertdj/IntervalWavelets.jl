@doc """
	weval(coef::Vector, res::Int) -> Vector

Evaluate the coefficients `coef` of the Haar scaling function basis on `[0,1]`.
The functions are evaluated at the dyadic rationals of resolution `res`.
"""->
function weval(coef::AbstractVector, res::Integer)
	@assert !isempty(coef)
	@assert res >= 0
	const Ncoef = length(coef)
	@assert ispow2(Ncoef)
	const J = Int( log2(Ncoef) )

	const mother_support = support("haar")
	const x = dyadic_rationals(mother_support, res)


	const Nx = length(x)
	y = zeros(Float64, Nx)
	for m in 1:Ncoef
		S = support(mother_support, J, m-1)
		for nx in 1:Nx
			if isinside(x[nx], S)
				@inbounds y[nx] += coef[m]*HaarScaling( x[nx], J, m-1 )
			end
		end
	end

	return x, y
end

@doc """
	weval(coef::Matrix, res::Int) -> Matrix

Evaluate the coefficients `coef` of the Haar scaling function basis on `[0,1]^2`.
The functions are evaluated at the dyadic rationals of resolution `res`.
"""->
function weval(coef::AbstractMatrix, res::Integer)
	@assert !isempty(coef)
	@assert (Ncoef = size(coef,1)) == size(coef,2)
	@assert res >= 0

	mother_support = support("haar")
	const x = dyadic_rationals(mother_support, res)
	const J = Int( log2(Ncoef) )
	const Nx = length(x)
	y = zeros(Float64, Nx, Nx)

	for my in 1:Ncoef, mx in 1:Ncoef
		Sx = support(mother_support, J, mx-1)
		Sy = support(mother_support, J, my-1)

		for ny in 1:Nx
			if !isinside(x[ny], Sy)
				continue
			end
			for nx in 1:Nx
				if isinside(x[nx], Sx)
					# TODO: Make Haar for 2D input?
					y[ny,nx] += coef[my,mx]*HaarScaling(x[nx], J, mx-1) * HaarScaling(x[ny], J, my-1)
				end
			end
		end
	end

	return y
end

@doc """
	weval(coef::Vector, N::Int, J::Int, res::Int)

Evaluate `coef` vector in the Daubechies `N` basis at scale `J` in the dyadic rationals of resolution `res`.

`J` must be so large that the support of the wavelet is contained in `[0,1]`.
"""->
function weval(coef::AbstractVector, N::Int, J::Integer, res::Integer)
	# TODO: J must be sufficiently large 

	x = dyadic_rationals(N, res)
	Nx = length(x)
	Ncoef = length(coef)
	# TODO: Compute J from coef?

	y = zeros(Float64, Nx)
	for n in 1:Nx
		for m in 1:Ncoef
			# TODO: Only include the functions that have x[nx] in their support
			@inbounds y[n] += coef[m]*HaarScaling( x[n], J, m-1 )
		end
	end

	return x, y
end


