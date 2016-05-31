module IntervalWavelets

import Wavelets: wavelet, WT

export
	# Types
	DaubSupport,
	InteriorFilter,
	BoundaryFilter,

	# Filter functions
	left,
	right,
	support,
	coef,
	van_moment,
	ifilter,
	bfilter,

	# Dyadic rationals
	dyadic_rationals,
	isuniform,
	x2index,
	index2x,
	x2index,
	index2x,

	# Scaling functions
	DaubScaling,

	# Reconstruction
	IntervalScaling,
	weval

const sqrt2 = sqrt(2)

include("Filters.jl")
include("Misc.jl")
include("DyadicPoints.jl")
include("Daubechies.jl")
include("BoundaryDaubechies.jl")
include("Quad.jl")
include("Reconstruction.jl")

end # module
