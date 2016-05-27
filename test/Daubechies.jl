using IntervalWavelets
using Base.Test

println("Testing interior scaling functions...")

#=
Orthogonality of integer translates of the internal scaling function
=#

EPS = 1e-3
R = 10
for N in 2:10
	C = ifilter(N)

	# Union of supports of all translated scaling functions
	S = support(C)
	extended_supp = DaubSupport( left(S), 2*right(S)-1 )
	x = dyadic_rationals( extended_supp, R )

	# Base scaling funtion
	y1 = zeros(x)
	phi = DaubScaling(C, R)
	@test_approx_eq_eps IntervalWavelets.l2norm(dyadic_rationals(S,R), phi) 1.0 EPS

	supp_index = [1:length(phi);]
	y1[supp_index] = phi

	# Integer translates
	y2 = similar(y1)
	for k in 1:right(S)-1
		fill!(y2, 0.0)
		y2[k*2^R-1+supp_index] = phi
		#@show IntervalWavelets.inner(x, y1, y2)
		@test_approx_eq_eps IntervalWavelets.inner(x, y1, y2) 0.0 EPS
	end
end

