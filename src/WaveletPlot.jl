module WaveletPlot

using Wavelets

export
	# Functions
	dyadic_rationals,
	x2index,
	index2x,
	support,
	DaubScaling

include("Types.jl")
include("DyadicPoints.jl")
include("Daubechies.jl")
#= include("BoundaryDaubechies.jl") =#

end # module
