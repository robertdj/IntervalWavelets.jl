using WaveletPlot
import WaveletPlot: inner, unit
using Base.Test

println("Testing reconstruction in wavelet basis...")

#=
- Scaling functions related to wavelets with p vanishing moments should reconstruct polynomials of degree p
- Unit vectors should reconstruct a single function
=#

R = 10
J = 4
S = R - J
DS = DaubSupport(0,1)

for p = 2:8
	for k = 1:2^J
		u = unit(2^J, k)
		xrecon, yrecon = weval(u, p, R)

		# The reconstructed function is dilated, scaled and maybe translated so
		# the non-zero part is extracted
		if 1 <= k <= p
			x, Y = DaubScaling(p, 'L', S)
			y = Y[:,k]

			sidx = 1
			eidx = x2index( 2.0^-J*(2*p-1), DS, R)
		elseif p < k <= 2^J-p
			x, y = DaubScaling(p, S, true)

			DS = DaubSupport(0,1)
			sidx = x2index( 2.0^(-J)*(-p+k), DS, R )
			eidx = x2index( 2.0^(-J)*(p+k-1), DS, R )
		elseif k > 2^J-p
			x, Y = DaubScaling(p, 'R', S)
			y = Y[:,2^J-k+1]

			sidx = x2index( 1-2.0^(-J)*(2*p-1), DS, R )
			eidx = x2index( 1, DS, R )
		end

		yrecon = yrecon[sidx:eidx]
		#@show k, norm(yrecon-sqrt(2.0^J)*y)
		@test_approx_eq yrecon sqrt(2.0^J)*y
	end
end

