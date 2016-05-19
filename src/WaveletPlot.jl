module WaveletPlot

import Wavelets: wavelet, WT

export
	# Types
	ScalingFilters,
	DaubSupport,
	BoundaryFilter,

	# Functions
	left,
	right,
	dyadic_rationals,
	isuniform,
	x2index,
	index2x,
	support,
	DaubScaling,
	waveparse,
	van_moment,
	ifilter,
	bfilter,
	weval

const sqrt2 = sqrt(2)

include("Filters.jl")
include("misc.jl")
include("DyadicPoints.jl")
include("Daubechies.jl")
include("BoundaryDaubechies.jl")
include("DilationTranslation.jl")
include("Quad.jl")
include("Reconstruction.jl")

end # module
