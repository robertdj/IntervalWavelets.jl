module IntervalWavelets

import OffsetArrays
import LinearAlgebra
import Wavelets

export 
    DyadicRational,
    DyadicRationalVector,
    Filter,
    InteriorFilter,

    coefficients,
    interior_filter,
    numerator,
    scale,
    support,
    vanishing_moments

const sqrt2 = sqrt(2)

include("DyadicRationals.jl")
include("Filters.jl")
include("Filters/Interior.jl")
include("Filters/Left.jl")
include("Filters/Right.jl")
include("InteriorScalingFunction.jl")
include("Misc.jl")

end # module
