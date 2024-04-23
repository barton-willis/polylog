# polylog

## Introduction

 The file `polylog.jl` has code for the numerical evaluation of the dilogarithm function on the entire complex plane. For a definition of this function, see the [Digital Library of Mathematical Functions](https://dlmf.nist.gov/25.12#E1)

The dilogarithm function has applications to statistical mechanics, Feynman integrals,
and algebraic number theory, to name a few.

 This Julia code uses a series representation from the article “The binomial transform of p-recursive sequences and the dilogarithm function,” [(arXiv)][def] by Stephanie Harshbarger and Barton Willis. This article was published by _Applications and Applied Mathematics_, in **vol.** 15, **Issue** 2 (December 2020), pp. 1025–1031.

## Features

- The same code is used for real and complex `binary16`, `binary32`, `binary64`, and `bigfloat` numbers.

- The method uses a recursively defined sequence that converges linearly on 
$\{z \in \{\mathbf{C} \mid \mathrm{Re}(z) \leq 1 \}$ to $\mathrm{Li}_2(z)$. The recursion for the sequence $\Phi$ is stable, meaning that asymptotically, all solutions to the recursion are subdominant to the solution that converges to the dilogarithm function. To gain a better convergence rate, this method uses various dilogarithm functional relations.

- Additionally, there is another method based on a sequence that converges for $\mathbf{C} \setminus [1,\infty)$ to $\mathrm{Li}_2(z)$. Asymptotically, the sequence $\Phi$ satisfies $\Phi_k \sim L + a \mu^k$, for $k \to \infty$, where the magnitude of $\mu$ is bounded by $1/\sqrt{3} \approx 0.577$ and $a \in \mathbf{C}$. Although this series has a somewhat better convergence rate than the series that converges on $\{z \in \mathbf{C} \mid \mathrm{Re}(z) \leq 1 \}$, this method is generally somewhat slower and the code is less tested.

- The focus of this code is _accuracy over speed_. To boost accuracy, the method uses Kahan summation. Also, it monitors the accuracy of the summation using a running error bound.

## Other

There is a standard Julia package [PolyLog.jl](https://juliapackages.com/p/polylog) for the numerical evaluation of polylogarithms of various orders, not just $\mathrm{Li}_2$. For binary64 numbers, the Julia package `PolyLog.jl` uses efficient rational function approximations, and its speed is _far_
better than `polylog.jl`. For example

~~~Julia
x = cis(pi/3);
@btime polylog2(x)
    4.529 μs (1 allocation: 32 bytes)
0.27415567780803785 + 1.0149416064096537im
@btime li2(x)
  206.140 ns (1 allocation: 32 bytes)
0.27415567780803807 + 1.0149416064096537im
~~~

For BigFloat numbers, the method in `polylog.jl` is _usually_ slightly faster than `PolyLog.jl`, but uses more memory. Here is one example

~~~Julia
setprecision(BigFloat,128);

@btime li2(convert(Complex{BigFloat},cis(pi/3)))
  4.622 ms (4940 allocations: 158.87 KiB)

0.2741556778080378663699490634254841513964 + 1.14941606409653637648270733876243611287im

@btime polylog2(convert(Complex{BigFloat},cis(pi/3)))

  2.444 ms (42803 allocations: 1.61 MiB)

0.2741556778080378663699490634254841514023 + 1.14941606409653637648270733876243611281im
~~~

 The file `tests.jl` has some tests of special values and some dilogarithm function identity based tests.

[def]: https://arxiv.org/pdf/1910.06928.pdf
