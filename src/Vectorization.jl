@vectorize_1arg Real HaarScaling

function HaarScaling(xi::AbstractArray, J::Int, k::Int)
	y = Array(Float64, size(xi))
	N = length(xi)
	for n = 1:N
		@inbounds y[n] = HaarScaling( xi[n], J, k )
	end
	return y
end



