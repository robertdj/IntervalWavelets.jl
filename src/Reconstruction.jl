@doc """
	weval(coef::Vector, res::Int) -> Vector

Evaluate the coefficients `coef` of the Haar scaling function basis on `[0,1]`.
The functions are evaluated at the dyadic rationals of resolution `res`.
"""->
function weval(coef::AbstractVector, res::Integer)
	Ncoef = length(coef)
	@assert ispow2(Ncoef)
	@assert (J = Int(log2(Ncoef)) ) <= res

	# Points in support
	recon_supp = DaubSupport(0,1)
	x = dyadic_rationals(recon_supp, res)

	# TODO: dilation -> Integer
	#= dilation = 2.0^(-J) =#
	dilation = 1/Ncoef
	sqrt_dil = sqrt(dilation)
	y = zeros(x)

	for k in 0:Ncoef-1
		# The indices of the x's in the support of phi_{J,k}
		start_idx = x2index( dilation*k, recon_supp, res )
		end_idx = x2index( dilation*(k+1), recon_supp, res )

		for nx in start_idx:end_idx
			phi = HaarScaling(x[nx]/dilation - k) / sqrt_dil
			@inbounds y[nx] += coef[k+1]*phi
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
	@assert ( J = Int(log2(Ncoef)) ) <= res

	recon_supp = DaubSupport(0,1)
	x = dyadic_rationals(recon_supp, res)

	Nx = length(x)
	y = zeros(Float64, Nx, Nx)
	dilation = 2.0^(-J)
	sqrt_dil = sqrt(dilation)

	for ky in 0:Ncoef-1, kx in 0:Ncoef-1
		startx_idx = x2index( dilation*kx, recon_supp, res )
		endx_idx = x2index( dilation*(kx+1), recon_supp, res )

		starty_idx = x2index( dilation*ky, recon_supp, res )
		endy_idx = x2index( dilation*(ky+1), recon_supp, res )

		for nx in startx_idx:endx_idx
			@inbounds phix = HaarScaling(x[nx]/dilation - kx) / sqrt_dil

			for ny in starty_idx:endy_idx
				@inbounds phiy = HaarScaling(x[ny]/dilation - ky) / sqrt_dil
				@inbounds y[ny,nx] += coef[ky+1,kx+1] * phix * phiy
			end
		end
	end

	return y
end

function weval(coef::AbstractVector, p::Integer, R::Integer)
	Ncoef = length(coef)

	# Points in support
	recon_supp = DaubSupport(0,1)
	x = dyadic_rationals(recon_supp, R)

	# Allocate output and a temporary array for the individual terms in
	# the reconstruction
	y = zeros(x)
	cur_phi = zeros(x)
	ZERO = zero( eltype(cur_phi) )

	# ----------------------------------------
	# Left side

	S = R - J
	lfilter = bfilter(p, 'L')
	int_filter = ifilter(p, true)
	lphi = DaubScaling(lfilter, int_filter, S)'
	sqrt_dil = sqrt(Ncoef)
	scale!(lphi, sqrt_dil)
	Nphi = size(lphi,1)

	coef_idx = 0
	for k = 0:p-1
		source_offset = Nphi*k + 1
		copy!( cur_phi, 1, lphi, source_offset, Nphi )
		BLAS.axpy!( coef[coef_idx+=1], cur_phi, y )
	end

	# ----------------------------------------
	# Interior

	iphi = DaubScaling(int_filter, S)
	scale!(iphi, sqrt_dil)

	inv_dil = 1/Ncoef
	for k in p:Ncoef-p-1
		dest_offset = x2index( inv_dil*(k-p+1), recon_supp, R )

		fill!( cur_phi, ZERO )
		copy!( cur_phi, dest_offset, iphi, 1, Nphi )
		BLAS.axpy!( coef[coef_idx+=1], cur_phi, y )
	end

	# ----------------------------------------
	# Right side

	rfilter = bfilter(p, 'R')
	rphi = DaubScaling(rfilter, int_filter, S)
	# rphi contains by in the "opposite" order: The first column is the
	# function closest to the boundary
	rphi = flipdim(rphi,1)'
	scale!(rphi, sqrt_dil)

	dest_offset = x2index( 1-inv_dil*(2*p-1), recon_supp, R )
	for k in 0:p-1
		source_offset = Nphi*k + 1
		# TODO: support/offset depends on k
		#= dest_offset = x2index( 1-inv_dil*(p+k), recon_supp, R ) =#

		fill!( cur_phi, ZERO )
		copy!( cur_phi, dest_offset, rphi, source_offset, Nphi )
		BLAS.axpy!( coef[coef_idx+=1], cur_phi, y )
	end

	return x, y
end

