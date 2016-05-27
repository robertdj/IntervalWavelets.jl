import IntervalWavelets: trapezquad
using Base.Test

# Test trapezquad
x = dyadic_rationals( DaubSupport(-1,1), 3 )
@test trapezquad( x, ones(x) ) == 2
@test trapezquad( x, x ) == 0

