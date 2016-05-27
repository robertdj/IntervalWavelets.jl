# WaveletPlot

[![Build Status](https://travis-ci.org/robertdj/WaveletPlot.jl.svg?branch=master)](https://travis-ci.org/robertdj/WaveletPlot.jl)
[![codecov.io](https://codecov.io/github/robertdj/WaveletPlot.jl/coverage.svg?branch=master)](https://codecov.io/github/robertdj/WaveletPlot.jl?branch=master)

*WaveletPlot* is a Julia package for computing ordinary/interior Daubechies scaling functions and the moment preserving boundary scaling functions of Cohen, Daubechies and Vial.
See the [enclosed document](doc/boundary_wavelets.pdf) for further description of these functions.


## Scaling functions

A Daubechies scaling function is defined by a filter that is related to the number of vanishing moments of the associated wavelet.
A filter for the ordinary/interior Daubechies scaling function "with" `p` vanishing moments is computed as

```julia
C = ifilter(p)
```

A second, boolean argument determines whether or not to use the linear phase (`true`) or minimum phase (`false`) scaling function; the default is `true`.

Except for the explicit [Haar wavelet](https://en.wikipedia.org/wiki/Haar_wavelet), Daubechies scaling functions can only be calculated at the [dyadic rationals](https://en.wikipedia.org/wiki/Dyadic_rational), i.e., points of the form k/2^R for R >= 0.
In *WaveletPlot* R is called the *resolution* of the dyadic rationals.
A Daubechies scaling function defined by a filter `C` is evaluated at all dyadic rationals of resolution R in its support with the command

```julia
y = DaubScaling(C, R)
```

Alternatively, use

```julia
x, y = DaubScaling(p, R)
```

that also returns the points where the function is evaluated.

For `p` vanishing moments there is also `p` left and `p` right scaling functions.
Their filters are available with as

```julia
L = bfilter(p, 'L')
R = bfilter(p, 'R')
```

If `B` denotes either `L` or `R` above and `C` is the interior filter, the `p` boundary functions evaluated at the dyadic rationals of resolution R are returned as columns in the matrix `Y` with the command

```julia
Y = DaubScaling(B, C, R)
```

Alternatively,

```julia
x, Y = DaubScaling(p, 'L', R)
x, Y = DaubScaling(p, 'R', R)
```

returns the points where the functions are evaluated, as for the interior scaling function.


## Compute representations

A representation of a 1D or 2D function in a scaling function basis is a vector or matrix of coefficients with respect to that basis.
The reconstruction of the function from this array is a linear combination of scaling functions.
In *WaveletPlot* reconstructions are computed on the unit interval.

Let `coeff` be an array with side(s) that is/are a power of 2.
To compute the reconstruction in the Haar basis in the dyadic rationals of resolution `R`, use

```julia
f = weval(coeff, R)
```

To compute the reconstruction in with Daubechies scaling functions "with" `p` vanishing moments at resolution `R`, use

```julia
f = weval(coeff, p, R)
```


## Miscellaneous

A lot of the testing is based on computing inner products in L^2([0,1]) and to do this a couple of functions are included:

- `trapezquad(x, y)`: Integrate `y` over `x` using the [trapezoidal rule](https://en.wikipedia.org/wiki/Trapezoidal_rule).
- `inner(x, y, z)`: The L2 inner product of `y` and `z` over `x`.
- `l2norm(x, y)`: The L2 norm of `y` over `x`.


## Installation

In Julia, simply run

```julia
Pkg.add("WaveletPlot")
```
