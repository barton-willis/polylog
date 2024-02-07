# polylog

## Introduction

 The file `polylog.jl` has code for the numerical evaluation of the dilogarithm function on the entire complex plane. For a definition of this function, see the [Digital Library of Mathematical Functions](https://dlmf.nist.gov/25.12#E1)

 The method is based on a series representation from the article "The binomial transform of p-recursive sequences and the dilogarithm function," [(arXiv)][def] by Stephanie Harshbarger and Barton Willis. This article was published by _Applications and Applied Mathematics_, in **vol.** 15, **Issue** 2 (December 2020), pp. 1025–1031.

## Features

- The same code is used for real and complex `binary16`, `binary32`, `binary64`, and `bigfloat` numbers.

- The method uses a sequence that converges linearly on the entire complex plane (except at one). Asymptotically, the approximating sequence $\Phi$ satisfies $\Phi_k \sim   L + a \mu^k$, for $k \to \infty$, where the magnitude of $\mu$ is bounded by $1/\sqrt{3} \approx 0.577$ and $a \in \mathbf{C}$. Further, the recursion for the sequence $\Phi$ is stable, meaning that asymptotically, all solutions to the recursion are
subdominant to the solution that converges
to the dilogarithm function.

- The focus of this code is _accuracy over speed_. To boost accuracy, the method uses Kahan summation. Also, it monitors the accuracy of the summation using a running error bound.

## Other

There is a standard Julia package [PolyLog.jl](https://juliapackages.com/p/polylog) for the numerical evaluation of polylogarithms. For binary64 numbers, the Julia package
`PolyLog.jl` uses efficient rational function approximations, and its speed is far
better than `polylog.jl`. For example

~~~Julia
x = cis(pi/3);
@btime polylog2(x)
  2.256 μs (21 allocations: 464 bytes)
  0.27415567780803785 + 1.0149416064096535im
@btime li2(x)
  139.277 ns (1 allocation: 32 bytes)
  0.27415567780803785 + 1.0149416064096537im
~~~

But for BigFloat numbers, the method in `polylog.jl` _is sometimes_ more efficient than `PolyLog.jl`. Here is one example

~~~Julia
setprecision(BigFloat,128);

@btime li2(convert(Complex{BigFloat},cis(pi/3)))
2.079 ms (4912 allocations: 156.58 KiB)
0.2741556778080378663699490634254841514023 + 1.14941606409653637648270733876243611287im

@btime polylog2(convert(Complex{BigFloat},cis(pi/3s)))
857.000 μs (22015 allocations: 871.88 KiB)
0.2741556778080378663699490634254841514082 + 1.14941606409653637648270733876243611281im
~~~

 The file `tests.jl` has some tests of special values. These tests make use of the standard Julia unit testing
format. The file `polylog.jl` has some identity-based tests, but these tests are not (yet) in the form of Julia unit testing.
So far, I have not attempted to build a Julia package.

[def]: https://arxiv.org/pdf/1910.06928.pdf
