module WaveletPlot

import Wavelets: wavelet, WT

const sqrt2 = sqrt(2)

export
	# Types
	ScalingFilters,
	BoundaryFilter,

	# Functions
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

include("Types.jl")
include("misc.jl")
include("DyadicPoints.jl")
include("Daubechies.jl")
include("BoundaryDaubechies.jl")
include("DilationTranslation.jl")
include("Filters.jl")
include("Quad.jl")
include("Reconstruction.jl")

end # module
