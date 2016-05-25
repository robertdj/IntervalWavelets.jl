#= using WaveletPlot =#
using Polynomials
#= using Base.Test =#

println("Testing reconstruction in wavelet basis...")

#=
- Scaling functions related to wavelets with p vanishing moments should reconstruct polynomials of degree p
=#

R = 10
J = 4
dilation = 2.0^J
sqrt_dil = sqrt(dilation)
coef = zeros(2^J)

for p = 2:4
	println(p)
	for m = 0:p+2
		# Project polynomial on scaling functions
		q = Poly( rand(m+1) )

		xL, yL = DaubScaling(p, 'L', R)
		for k = 1:p
			qeval = q( xL/dilation )
			coef[k] = inner(xL, yL[:,k], qeval) / sqrt_dil
		end

		xI, yI = DaubScaling(p, R, true)
		for k = p:2^J-p-1
			qeval = q( (xI+k)/dilation )
			coef[k+1] = inner(xI, yI, qeval) / sqrt_dil
		end

		xR, yR = DaubScaling(p, 'R', R)
		qeval = q( 1 + xR/dilation )
		for k = 1:p
			coef[end+1-k] = inner(xR, yR[:,k], qeval) / sqrt_dil
		end

		# Reconstruct polynomial from coefficients
		x, y = weval(coef, p, R)
		qeval = q(x)
		@show norm(qeval - y)
	end
end

