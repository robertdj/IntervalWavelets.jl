module IntervalWavelets

using OffsetArrays
using Compat
using RecipesBase

import Compat: view, String
import Wavelets: wavelet, WT

export
	# Types
	DaubSupport,
	Interval,
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
	DyadicRationalsVector,
	dyadic_rationals,
	isuniform,
	x2index,
	index2x,
	x2index,
	index2x,

	# Scaling functions
	HaarScaling,
	DaubScaling,

	# Reconstruction
	IntervalScaling,
	weval

const sqrt2 = sqrt(2)

include("DyadicPoints.jl")
include("Filters.jl")
include("Misc.jl")
include("Daubechies.jl")
include("Vectorization.jl")
include("BoundaryDaubechies.jl")
include("Quad.jl")
include("Reconstruction.jl")

end # module
