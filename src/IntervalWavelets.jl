module IntervalWavelets

import OffsetArrays
import LinearAlgebra
using RecipesBase
import Wavelets

export 
    DyadicRational,
    DyadicRationalVector,
    Filter,
    InteriorFilter,

    coefficients,
    interior_filter,
    interior_scaling_function,
    numerator,
    reduce,
    scale,
    support,
    support_boundaries,
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
