module IntervalWavelets

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
	coef,
	DaubScaling,
	waveparse,
	van_moment,
	ifilter,
	bfilter,
	weval,
	x2index,
	index2x

const sqrt2 = sqrt(2)

include("Filters.jl")
include("Misc.jl")
include("DyadicPoints.jl")
include("Daubechies.jl")
include("BoundaryDaubechies.jl")
include("Quad.jl")
include("Reconstruction.jl")

end # module
