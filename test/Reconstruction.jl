using IntervalWavelets
import IntervalWavelets: inner, unit
using Base.Test

println("Testing reconstruction in wavelet basis...")

#=
- Scaling functions related to wavelets with p vanishing moments should reconstruct polynomials of degree p
- Unit vectors should reconstruct a single function
=#

R = 7
J = 4
S = R - J
DS = DaubSupport(0,1)


# ------------------------------------------------------------------------
# 1D

for p = 2:8
	for k = 0:2^J-1
		u = unit(2^J, k+1)
		xrecon, yrecon = weval(u, p, R, DaubSupport(0,1))

		y = IntervalScaling(p, k, J, R)
		#@show k, norm(yrecon - y)
		@test_approx_eq yrecon y
	end
end


# ------------------------------------------------------------------------
# 2D

for p = 2:8
	p = 2
	# The tests pass for for all l and k, but takes forever w/o inlining
	k = rand( [0:2^J-1;] )
	l = rand( [0:2^J-1;] )

	u1 = unit(2^J, k+1)
	u2 = unit(2^J, l+1)
	U = u1*u2'
	yrecon = weval(U, p, R, DaubSupport(0,1))

	y1 = IntervalScaling(p, k, J, R)
	y2 = IntervalScaling(p, l, J, R)
	Y = y1*y2'
	#@show k, l, norm(yrecon - Y)
	@test_approx_eq yrecon Y
end

