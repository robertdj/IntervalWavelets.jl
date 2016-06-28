@doc """
	weval(coef, wavename, res)

Evaluate the coefficients `coef` of the `wavename` scaling function basis on `[0,1]` or `[0,1]^2` (depending on whether `coef` is a vector or a matrix).

The functions are evaluated at the dyadic rationals of resolution `res`.
"""->
function weval(coef, wavename::AbstractString, res::Integer)
	vm = van_moment(wavename)
	if vm == 1
		return weval(coef, res)
	elseif vm >= 2
		return weval(coef, vm, res)
	else
		error()
	end
end

@doc """
	weval(coef::Vector, res::Int, supp) -> x, y

Evaluate the coefficients `coef` of the Haar scaling function basis on the DaubSupport interval `supp`.
The functions are evaluated at the dyadic rationals of resolution `res`.
"""->
function weval(coef::AbstractVector, res::Integer, supp::DaubSupport)
	Ncoef = length(coef)
	lsupp = length(supp)
	Jpow = div(Ncoef, lsupp)
	Ncoef == lsupp * Jpow || throw(AssertionError())
	ispow2(Jpow) || throw(AssertionError())
	(J = Int(log2(Jpow)) ) <= res || throw(AssertionError())

	# Points in support
	x = dyadic_rationals(supp, res)

	dilation = 1/Jpow
	sqrt_dil = sqrt(dilation)
	y = zeros(x)

	translations = L*Jpow + (0:Ncoef-1)
	coef_idx = 0
	for k in translations
		# The indices of the x's in the support of phi_{J,k}
		start_idx = x2index( dilation*k, supp, res )
		end_idx = x2index( dilation*(k+1), supp, res )

		coef_idx+=1
		@inbounds for nx in start_idx:end_idx
			phi = HaarScaling(x[nx]/dilation - k) / sqrt_dil
			y[nx] += coef[coef_idx]*phi
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
	res >= 0 || throw(DomainError())
	(Ncoef = size(coef,1)) == size(coef,2) || throw(DimensionMismatch())
	ispow2(Ncoef) || throw(AssertionError())
	(J = Int(log2(Ncoef)) ) <= res || throw(DomainError())

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
	ispow2(Ncoef) || throw(AssertionError())
	2 <= p <= 8 || throw(AssertionError())
	log2(2*p-1) <= ( J = Int(log2(Ncoef)) ) <= R || throw(AssertionError())

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
	(Ncoef = size(coef,1)) == size(coef,2) || throw(DimensionMismatch())
	ispow2(Ncoef) || throw(AssertionError())
	2 <= p <= 8 || throw(AssertionError())
	log2(2*p-1) <= ( J = Int(log2(Ncoef)) ) <= R || throw(AssertionError())

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

	if p <= m < twopowJ - p
		unsafe_copy!( phi, 1, Φ, sz[1]*p+1, sz[1] )
	elseif m < p
		unsafe_copy!( phi, 1, Φ, sz[1]*m+1, sz[1] )
	else
		unsafe_copy!( phi, 1, Φ, sz[1]*(sz[2]-twopowJ+m)+1, sz[1] )
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
	1 <= p <= 8 || throw(AssertionError())
	R >= 0 || throw(DomainError())

	lfilter = bfilter(p, 'L')
	int_filter = ifilter(p, true)
	lphi = DaubScaling(lfilter, int_filter, R)'

	iphi = DaubScaling(int_filter, R)

	rfilter = bfilter(p, 'R')
	rphi = DaubScaling(rfilter, int_filter, R)'

	hcat(lphi, iphi, flipdim(rphi,2))
end

@doc """
	IntervalScaling(p::Int, k::Int, J::Int, R::Int) -> x, y

On the interval [0,1] there are 2^J functions at scale `J` ordered such
that the first `p` are left boundary scaling functions, the last `p` are
right boundary scaling functions and the middle 2^J-2p are translated
interior scaling functions.

`IntervalScaling` returns the `k`'th function computed at the dyadic rationals of resolution `R`.

Note that the right boundary scaling functions are ordered differently
than the output of `DaubScaling`:
The output of 

	IntervalScaling(p, 2^J, J, R)

is (a translated, dilated and scaled version of) the first column of `Y`
in

	x, Y = DaubScaling(p, 'R', R)
"""->
function IntervalScaling(p::Integer, k::Integer, J::Integer, R::Integer)
	2 <= p <= 8 || throw(AssertionError())
	log2(2*p-1) <= J <= R || throw(AssertionError("The scale is too small compared to p or too large compared to R"))
	0 <= k < (Ny = 2^J) || throw(AssertionError("This interval wavelet does not exist"))

	IF = ifilter(p)
	RmJ = R - J
	DS = DaubSupport(0,1)
	if 0 <= k < p
		BF = bfilter(p,'L')
		Y = DaubScaling(BF, IF, RmJ)
		y = vec( Y[k+1,:] )

		# The indices of the support
		start_idx = x2index(0, DS, R)
		end_idx = x2index( 2.0^-J*(2*p-1), DS, R)
	elseif p <= k < Ny-p
		y = DaubScaling(IF, RmJ)

		start_idx = x2index( 2.0^-J*(-p+k+1), DS, R )
		end_idx = x2index( 2.0^-J*(p+k), DS, R )
	else
		BF = bfilter(p,'R')
		Y = DaubScaling(BF, IF, RmJ)
		y = vec( Y[Ny-k,:] )

		start_idx = x2index( 1-2.0^-J*(2*p-1), DS, R )
		end_idx = x2index( 1, DS, R )
	end

	phi = zeros(Float64, 2^R+1)
	scale!(y, sqrt(Ny))
	phi[start_idx:end_idx] = y

	return phi
end

@doc """
	IntervalScaling(p::Int, J::Int, R::Int) -> matrix

Returns all 2^J scaling functions of scale `J` at resolution `R` on
[0,1] as a `2^R+1`-by-2^J`.

**Note:** 
You should probably only use this function for small values of `J`.
"""->
function IntervalScaling(p::Integer, J::Integer, R::Integer)

	Y = Array{Float64}(2^R+1, 2^J)
	for k in 0:2^J-1
		Y[:,k+1] = IntervalScaling(p, k, J, R)
	end

	return Y
end

