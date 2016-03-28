using WaveletPlot
using Wavelets
using Base.Test

# ------------------------------------------------------------
# Translation between index and values for dyadic rationals

vm = 2
C = wavelet( WT.Daubechies{vm}() ).qmf
#= F = scalingfilters(vm) =#
#= IS, BS = support(F) =#
IS = support(C)

#= for support_type in [:IS, :BS] =#
for support_type in [:IS]
	@eval begin
		# Level 0 = integers
		supportx = dyadic_rationals( $support_type, 0 )
		for x in round(Int, supportx)
			@test supportx[ x2index(x,$support_type) ] == x
		end

		for idx = 1:length(supportx)
			@test supportx[idx] == index2x(idx,$support_type)
		end

		# Higher level
		R = 2
		supportx = dyadic_rationals( $support_type, R )
		for x in supportx
			@test supportx[ x2index(x,$support_type,R) ] == x
		end

		for idx = 1:length(supportx)
			@test supportx[idx] == index2x(idx,$support_type,R)
		end
	end
end


