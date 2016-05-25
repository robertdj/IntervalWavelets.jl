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

	dilation = 1/Ncoef
	sqrt_dil = sqrt(dilation)
	y = zeros(x)

	for k in 0:Ncoef-1
		# The indices of the x's in the support of phi_{J,k}
		start_idx = x2index( dilation*k, recon_supp, res )
		end_idx = x2index( dilation*(k+1), recon_supp, res )

		@inbounds for nx in start_idx:end_idx
			phi = HaarScaling(x[nx]/dilation - k) / sqrt_dil
			y[nx] += coef[k+1]*phi
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

	for ky in 0:Ncoef-1
		starty_idx = x2index( dilation*ky, recon_supp, res )
		endy_idx = x2index( dilation*(ky+1), recon_supp, res )

		for kx in 0:Ncoef-1
			startx_idx = x2index( dilation*kx, recon_supp, res )
			endx_idx = x2index( dilation*(kx+1), recon_supp, res )

			for nx in startx_idx:endx_idx
				@inbounds phix = HaarScaling(x[nx]/dilation - kx) / sqrt_dil

				@inbounds for ny in starty_idx:endy_idx
					phiy = HaarScaling(x[ny]/dilation - ky) / sqrt_dil
					y[ny,nx] += coef[ky+1,kx+1] * phix * phiy
				end
			end
		end
	end

	return y
end

function weval(coef::AbstractVector, p::Integer, R::Integer)
	Ncoef = length(coef)
	@assert ispow2(Ncoef)
	@assert 2 <= p <= 8
	@assert log2(2*p-1) <= ( J = Int(log2(Ncoef)) ) <= R

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

	lphi = DaubScaling(lfilter, int_filter, S)
	Nphi = size(lphi,2)

	sqrt_dil = sqrt(Ncoef)
	scale!(lphi, sqrt_dil)

	coef_idx = 0
	for k = 0:p-1
		unsafe_copy!( cur_phi, 1, 1, lphi, k+1, p, Nphi )
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
		unsafe_copy!( cur_phi, dest_offset, iphi, 1, Nphi )
		BLAS.axpy!( coef[coef_idx+=1], cur_phi, y )
	end

	# ----------------------------------------
	# Right side

	rfilter = bfilter(p, 'R')
	rphi = DaubScaling(rfilter, int_filter, S)
	scale!(rphi, sqrt_dil)

	dest_offset = x2index( 1-inv_dil*(2*p-1), recon_supp, R )
	# rphi is in "opposite" order, i.e., the first row is the 
	# function closest to the boundary
	coef_idx = Ncoef + 1
	for k in 0:p-1
		# TODO: support/offset depends on k
		#= dest_offset = x2index( 1-inv_dil*(p+k), recon_supp, R ) =#

		fill!( cur_phi, ZERO )
		unsafe_copy!( cur_phi, dest_offset, 1, rphi, k+1, p, Nphi )
		BLAS.axpy!( coef[coef_idx-=1], cur_phi, y )
	end

	return x, y
end


function weval2(coef::AbstractVector, p::Integer, R::Integer)
	Ncoef = length(coef)
	@assert ispow2(Ncoef)
	@assert 2 <= p <= 8
	@assert log2(2*p-1) <= ( J = Int(log2(Ncoef)) ) <= R

	# Collect mother scaling functions
	S = R - J
	lfilter = bfilter(p, 'L')
	int_filter = ifilter(p, true)
	lphi = DaubScaling(lfilter, int_filter, S)'

	iphi = DaubScaling(int_filter, S)

	rfilter = bfilter(p, 'R')
	rphi = DaubScaling(rfilter, int_filter, S)'

	PHI = hcat(lphi, iphi, flipdim(rphi,2))
	sqrt_dil = sqrt(Ncoef)
	scale!(PHI, sqrt_dil)

	# Allocate output and temporary array
	y = zeros(2^R + 1)
	Nphi = length(iphi)
	phi = Array{Float64}(Nphi)

	# Translation with 1 is a fixed number of indices
	kinc = x2index( 1/Ncoef, DaubSupport(0,1), R ) - x2index( 0, DaubSupport(0,1), R )
	slice_idx = 1:Nphi

	sz = size(PHI)
	for k in 0:Ncoef-1
		unsafe_DaubScaling!(phi, k, PHI, Ncoef, sz, p)

		if p <= k <= Ncoef - p
			slice_idx += kinc
		end
		z = slice(y, slice_idx)
		BLAS.axpy!( coef[k+1], phi, z )
	end

	recon_supp = DaubSupport(0,1)
	x = dyadic_rationals(recon_supp, R)
	return x, y
end

function weval(coef::AbstractMatrix, p::Integer, R::Integer)
	@assert (Ncoef = size(coef,1)) == size(coef,2)
	@assert ispow2(Ncoef)
	@assert 2 <= p <= 8
	@assert log2(2*p-1) <= ( J = Int(log2(Ncoef)) ) <= R

	# Collect mother scaling functions
	S = R - J
	lfilter = bfilter(p, 'L')
	int_filter = ifilter(p, true)
	lphi = DaubScaling(lfilter, int_filter, S)'

	iphi = DaubScaling(int_filter, S)

	rfilter = bfilter(p, 'R')
	rphi = DaubScaling(rfilter, int_filter, S)'

	PHI = hcat(lphi, iphi, flipdim(rphi,2))
	sqrt_dil = sqrt(Ncoef)
	scale!(PHI, sqrt_dil)

	# Allocate output and temporary arrays 
	y = zeros(2^R+1, 2^R+1)
	Nphi = length(iphi)
	phi = Array{Float64}(Nphi, Nphi)
	phi1 = Array{Float64}(Nphi)
	phi2 = similar(phi1)

	# Translation is a fixed number of indices
	kinc = x2index( 1/Ncoef, DaubSupport(0,1), R ) - x2index( 0, DaubSupport(0,1), R )

	slicex_idx = 1:Nphi
	sz = size(PHI)
	for kx in 0:Ncoef-1
		if p <= kx <= Ncoef - p
			slicex_idx += kinc
		end

		slicey_idx = 1:Nphi
		for ky in 0:Ncoef-1
			# TODO: Manually?
			unsafe_DaubScaling!(phi, phi1, phi2, (ky,kx), PHI, Ncoef, sz, p)

			if p <= ky <= Ncoef - p
				slicey_idx += kinc
			end

			z = slice(y, slicey_idx, slicex_idx)
			BLAS.axpy!( coef[ky+1,kx+1], phi, z )
		end
	end

	return y
end

@doc """
Overwrite `phi` with the `m`'th scaling function from a matrix `Φ` where 
the first `p` columns are the left scaling functions, the last `p` columns 
are the right scaling functions and the middle column is the interior
scaling function.
"""->
function unsafe_DaubScaling!( phi::DenseVector, m::Integer, Φ::DenseMatrix, 
	twopowJ::Integer, sz::Tuple=size(Φ), p::Integer=div(sz[2]-1,2) )

	# TODO: Is it better/faster if Φ is transposed?
	if m < p
		unsafe_copy!( phi, 1, Φ, sz[1]*m+1, sz[1] )
	elseif m >= twopowJ - p
		unsafe_copy!( phi, 1, Φ, sz[1]*(sz[2]-twopowJ+m)+1, sz[1] )
	else
		unsafe_copy!( phi, 1, Φ, sz[1]*p+1, sz[1] )
	end
end

@doc """
Overwrite `phi` with the outer product of `phi1` and `phi1` after they
have been overwritten with the `m[1]`'th and `m[2]`'th column of `Φ`, respectively.
"""->
function unsafe_DaubScaling!( phi::DenseMatrix, phi1::DenseVector, phi2::DenseVector,
	m::Tuple, Φ::DenseMatrix, twopowJ::Integer, sz::Tuple=size(Φ), p::Integer=div(sz[2]-1,2) )

	unsafe_DaubScaling!(phi1, m[1], Φ, twopowJ, sz, p)
	unsafe_DaubScaling!(phi2, m[2], Φ, twopowJ, sz, p)

	for j in 1:sz[1], i in 1:sz[1]
		@inbounds phi[i,j] = phi1[i] * phi2[j]
	end
end

