function HaarScaling(xi::AbstractArray, J::Int, k::Int)
	y = Array{Float64}(size(xi))
	for idx in eachindex(xi)
		@inbounds y[idx] = HaarScaling( xi[idx], J, k )
	end
	return y
end



