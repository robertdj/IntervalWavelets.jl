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


# ------------------------------------------------------------------------

function weval(coef::AbstractVector, p::Integer, R::Integer)
	Ncoef = length(coef)
	@assert ispow2(Ncoef)
	@assert 2 <= p <= 8
	@assert log2(2*p-1) <= ( J = Int(log2(Ncoef)) ) <= R

	# Collect mother scaling functions
	PHI = allDaubScaling(p, R-J)
	sqrt_dil = sqrt(Ncoef)
	scale!(PHI, sqrt_dil)

	# Allocate output and temporary array
	y = zeros(2^R + 1)
	Nphi = size(PHI,1)
	phi = Array{Float64}(Nphi)

	# Translation with 1 is a fixed number of indices
	kinc = x2index( 1/Ncoef, DaubSupport(0,1), R ) - 1
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
	PHI = allDaubScaling(p, R-J)
	sqrt_dil = sqrt(Ncoef)
	scale!(PHI, sqrt_dil)

	# Allocate output and temporary arrays 
	y = zeros(2^R+1, 2^R+1)
	Nphi = size(PHI,1)
	phi = Array{Float64}(Nphi, Nphi)
	phix = Array{Float64}(Nphi)
	phiy = similar(phix)

	# Translation is a fixed number of indices
	kinc = x2index( 1/Ncoef, DaubSupport(0,1), R ) - 1

	slicex_idx = 1:Nphi
	sz = size(PHI)
	for kx in 0:Ncoef-1
		if p <= kx <= Ncoef - p
			slicex_idx += kinc
		end

		unsafe_DaubScaling!(phix, kx, PHI, Ncoef, sz, p)
		slicey_idx = 1:Nphi
		for ky in 0:Ncoef-1
			unsafe_DaubScaling!(phiy, ky, PHI, Ncoef, sz, p)
			for j in 1:sz[1], i in 1:sz[1]
				@inbounds phi[i,j] = phiy[i] * phix[j]
			end

			if p <= ky <= Ncoef - p
				slicey_idx += kinc
			end

			z = slice(y, slicey_idx, slicex_idx)
			# BLAS.axpy! is using *a lot* of memory
			for j in 1:Nphi, i in 1:Nphi
				@inbounds z[i,j] += coef[ky+1,kx+1] * phi[i,j]
			end
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

	# TODO: Let interior be the first
	if m < p
		unsafe_copy!( phi, 1, Φ, sz[1]*m+1, sz[1] )
	elseif m >= twopowJ - p
		unsafe_copy!( phi, 1, Φ, sz[1]*(sz[2]-twopowJ+m)+1, sz[1] )
	else
		unsafe_copy!( phi, 1, Φ, sz[1]*p+1, sz[1] )
	end
end

@doc """
	allDaubScaling(p::Integer, R::Integer)

Collect the `p` left scaling functions, the interior scaling function
and the `p` right scaling functions as columns in a matrix (in that
order).

All functions are computed at resolution `R`.
"""->
function allDaubScaling(p::Integer, R::Integer)
	@assert 0 <= p <= 8
	@assert R >= 0

	lfilter = bfilter(p, 'L')
	int_filter = ifilter(p, true)
	lphi = DaubScaling(lfilter, int_filter, R)'

	iphi = DaubScaling(int_filter, R)

	rfilter = bfilter(p, 'R')
	rphi = DaubScaling(rfilter, int_filter, R)'

	hcat(lphi, iphi, flipdim(rphi,2))
end

