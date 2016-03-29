# ------------------------------------------------------------
# Dilation and translation

function HaarScaling(xi::Real, J::Int)
	2.0^(J/2)*HaarScaling(2^J*xi)
end

function HaarScaling(xi::Real, J::Int, k::Int)
	@assert 0 <= k < (dil_fact = 2.0^J) "Translation is too large"
	2.0^(J/2)*HaarScaling(dil_fact*xi - k)
end

