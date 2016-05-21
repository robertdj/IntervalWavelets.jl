@doc """
	weval(coef::Vector, res::Int) -> Vector

Evaluate the coefficients `coef` of the Haar scaling function basis on `[0,1]`.
The functions are evaluated at the dyadic rationals of resolution `res`.
"""->
function weval(coef::AbstractVector, res::Integer)
	@assert res >= 0
	Ncoef = length(coef)
	@assert ispow2(Ncoef)
	J = Int( log2(Ncoef) )

	# Points in support
	recon_supp = DaubSupport(0,1)
	x = dyadic_rationals(recon_supp, res)

	dilation = 2.0^(-J)
	y = zeros(x)

	for k in 0:Ncoef-1
		# The indices of the x's in the support of phi_{J,k}
		start_idx = x2index( dilation*k, recon_supp, res )
		end_idx = x2index( dilation*(k+1), recon_supp, res )

		for nx in start_idx:end_idx
			@inbounds y[nx] += coef[k+1]*HaarScaling( x[nx], J, k )
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
	@assert res >= 0
	@assert (Ncoef = size(coef,1)) == size(coef,2)
	@assert ispow2(Ncoef)
	@assert (J = Int(log2(Ncoef)) ) <= res

	recon_supp = DaubSupport(0,1)
	x = dyadic_rationals(recon_supp, res)

	Nx = length(x)
	dilation = 2.0^(-J)
	y = zeros(Float64, Nx, Nx)

	for ky in 0:Ncoef-1, kx in 0:Ncoef-1
		startx_idx = x2index( dilation*kx, recon_supp, res )
		endx_idx = x2index( dilation*(kx+1), recon_supp, res )

		starty_idx = x2index( dilation*ky, recon_supp, res )
		endy_idx = x2index( dilation*(ky+1), recon_supp, res )

		for ny in starty_idx:endy_idx
			for nx in startx_idx:endx_idx
				phix = HaarScaling(x[nx]/dilation - kx) / sqrt(dilation)
				phiy = HaarScaling(x[ny]/dilation - ky) / sqrt(dilation)
				@inbounds y[ny,nx] += coef[ky+1,kx+1] * phix * phiy
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


