using IntervalWavelets
using Base.Test

println("Testing boundary scaling functions...")

#=
The boundary wavelets are mutually orthonormal and orthogonal to the interior wavelets.
=#

EPS = 1e-3
R = 11

# ------------------------------------------------------------
# Boundary to boundary

for side in ['L'; 'R'], p in 2:8
	#@show side, p
	x, Y = DaubScaling(p, side, R)

	F = bfilter(p, side)
	S = support(F)

	for k in 1:p
		#@show IntervalWavelets.l2norm(x, Y[:,k])
		@test_approx_eq_eps IntervalWavelets.l2norm(x, Y[:,k]) 1.0 10*EPS

		for l in k+1:p
			intersect_supp = intersect( support(F,k-1), support(F,l-1) )
			xx = dyadic_rationals(intersect_supp, R)
			startx = x2index( left(intersect_supp), S, R )
			stopx = x2index( right(intersect_supp), S, R )

			#@show IntervalWavelets.inner(xx, Y[startx:stopx,k], Y[startx:stopx,l])
			@test_approx_eq_eps IntervalWavelets.inner(xx, Y[startx:stopx,k], Y[startx:stopx,l]) 0.0 EPS
		end
	end
end


# ------------------------------------------------------------
# Boundary to interior
# The interior function are translated such that there is a non-empty
# overlap with support of the boundary function.

EPS = 1e-5

for side in ['L'; 'R'], p in 2:8
	C = ifilter(p,true)
	interior = DaubScaling(C,R)
	suppC = support(C)
	int_index = [0:length(interior)-1;]

	B = bfilter(p, side)
	bound = DaubScaling(B, C, R)'

	# Boundary functions
	for k in 0:p-1
		# Translations of interior functions
		suppB = support(B,k)
		translation_index = (side == 'L' ? (p:2*p+k-2) : -2*p-k+1:-p-1)

		for l in translation_index
			suppI = intersect( suppC+l, suppB )
			x = dyadic_rationals(suppI, R)

			startB = x2index( left(suppI), suppB, R )
			stopB = x2index( right(suppI), suppB, R )

			startC = x2index( left(suppI), suppC+l, R )
			stopC = x2index( right(suppI), suppC+l, R )

			#@show IntervalWavelets.inner(x, interior[startC:stopC], bound[startB:stopB,k+1])
			@test_approx_eq_eps IntervalWavelets.inner(x, interior[startC:stopC], bound[startB:stopB,k+1]) 0.0 EPS
		end
	end
end

