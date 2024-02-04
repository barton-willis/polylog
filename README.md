# polylog

## Introduction

 The file `polylog.jl` has code for the numerical evaluation of the dilogarithm function on the entire complex plane. For a definition of this function, see the [Digital Library of Mathematical Functions](https://dlmf.nist.gov/25.12#E1)

 The method is based on a series representation from the article "The binomial transform of p-recursive sequences and the dilogarithm function," [(arxiv)][def] by Stephanie Harshbarger and Barton Willis. This article was published by _Applications and Applied Mathematics_, in **vol.** 15, **Issue** 2 (December 2020), pp. 1025–1031.

## Features

- The same code is used for real and complex `binary16`, `binary32`, `binary64`, and `bigfloat` numbers.

- The method uses a sequence that converges linearly on the entire complex plane (except at one). Specifically the sequence $\Phi$ satisfies $\Phi_k \sim   L + \mu^k$,
for $k \to \infty$, where the magnitude of $\mu$ is bounded by $1/\sqrt{3} \approx 0.577$.

- The focus of this code is _accuracy over speed_. To boost accuracy, the method uses Kahan summation. Also, it monitors the accuracy of the summation using a running error bound.

## Other

- There is a standard Julia package [PolyLog.jl](https://juliapackages.com/p/polylog) for the numerical evaluation of polylogarithms, but I have not compared the methods.

- The file `tests.jl` has some tests of special values. These tests make use of the standard Julia unit testing
format. The file `polylog.jl` has some identity-based tests, but these tests are not (yet) in the form of Julia unit testing.

- So far, I have not attempted to build a Julia package.

[def]: https://arxiv.org/pdf/1910.06928.pdf
