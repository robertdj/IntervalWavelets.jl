module IntervalWavelets

import OffsetArrays
import LinearAlgebra
import Wavelets

using RecipesBase

export 
    DyadicRational,
    DyadicRationalVector,
    Filter,
    InteriorFilter,
    BoundaryFilter,
    Sides,
    LEFT,
    RIGHT,

    boundary_filters,
    boundary_scaling_functions,
    coefficients,
    filters,
    interior_filter,
    interior_scaling_function,
    numerator,
    reduce,
    scale,
    side,
    support,
    support_boundaries,
    vanishing_moments

const sqrt2 = sqrt(2)

include("DyadicRationals.jl")
include("Enums.jl")
include("Filters.jl")
include("Filters/Interior.jl")
include("Filters/Left.jl")
include("Filters/Right.jl")
include("InteriorScalingFunction.jl")
include("BoundaryScalingFunction.jl")
include("Misc.jl")

end # module
