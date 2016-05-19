using WaveletPlot
using Base.Test

println("Testing boundary scaling functions...")

#=
The boundary wavelets are mutually orthonormal and orthogonal to the interior wavelets.
=#

EPS = 1e-3
R = 11

# ------------------------------------------------------------
# Boundary to boundary

for side in ['L'], p in 2:8
	x, Y = DaubScaling(p, side, R)

	for k in 1:p
		#@show WaveletPlot.l2norm(x, Y[:,k])
		@test_approx_eq_eps WaveletPlot.l2norm(x, Y[:,k]) 1.0 10*EPS
		for l in k+1:p
			#@show WaveletPlot.inner(x, Y[:,k], Y[:,l])
			@test_approx_eq_eps WaveletPlot.inner(x, Y[:,k], Y[:,l]) 0.0 EPS
		end
	end
end


# ------------------------------------------------------------
# Boundary to interior
# Interior function are translated such that the left support endpoint is at
# minimum 1 and until there's no overlap between the interior support
# and boundary support

EPS = 1e-5

for side in ['L'], p in 2:8
	C = ifilter(p,true)
	interior = DaubScaling(C,R)
	int_index = [0:length(interior)-1;]

	bound = DaubScaling(p, side, R)[2]

	# Boundary functions
	for k in 0:p-1
		# Translations of interior functions
		for l in 1:p+k-1
			supp_union = DaubSupport(0, 3*p-2+k)
			x = dyadic_rationals(supp_union, R)

			int_padded = zeros(x)
			int_supp = x2index(l, supp_union, R) + int_index
			int_padded[int_supp] = interior

			bound_padded = zeros(x)
			bound_padded[1:size(bound,1)] = bound[:,k+1]

			#@show WaveletPlot.inner(x, int_padded, bound_padded)
			@test_approx_eq_eps WaveletPlot.inner(x, int_padded, bound_padded) 0.0 EPS
		end
	end
end

