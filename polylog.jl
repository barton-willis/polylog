# Barton Willis, 2019, 2024

# This work is licensed under the CC0 1.0 Universal license.

# Julia code for the evaluation of the dilogarithm. The method is
# based on "The binomial transform of p-recursive sequences and 
# dilogarithm function," by Stephanie Harshbarger and Barton Willis.
# https://arxiv.org/pdf/1910.06928.pdf

# Extend eps to complex numbers.
import Base.eps
import Base.BigFloat
eps(::Type{Complex{T}}) where T <: AbstractFloat = eps(T)

import Base.precision
precision(::Type{Complex{T}}) where T <: AbstractFloat = precision(T)
# Is there a better way to do this?  Something like how pi is defined?
"""
  zeta2(T::Type)

Return the value of `pi^2/6` rounded to the type `T`

Examples:
```
julia> zeta2(Float64)
1.6449340668482264
julia> zeta2(Complex{Float32})
1.644934f0 + 0.0f0im
```
"""
function zeta2(T::Type)
    if T==Float64 || T==Complex{Float64}
        reinterpret(Float64, 4610086943623153619)
    elseif T==Float32 || T==Complex{Float32}
        reinterpret(Float32, convert(Int32, 1070763315))
    elseif T==Float16 || T==Complex{Float16}
        reinterpret(Float16, convert(Int16, 16020))
    else 
        convert(T,pi)*(convert(T,pi)/convert(T,6))
    end
end

"""
A complex natural logarithm function.

  clog(x::Number)

When x is in (0,infinity) return log(x); otherwise, convert x to 
a complex number and then dispatch log on x. This function is _not_
intended to be a user-level function.
"""
function clog(x::Number)
    # careful: for example 0.6 + 0.0im tests as real, but we don't want to send
    # 0.6 + 0.0im to log.
    if iszero(imag(x)) && real(x) > 0
        log(real(x))
    elseif iszero(imag(x))
        log(real(-x)) + convert(typeof(x),pi)*im
    else
        log(Complex(x))
    end
end

function clog1p(x)
    if iszero(imag(x)) && real(x) > 0
        log1p(x)
    else 
        log1p(Complex(x))
    end
end

"""
    KahanSum(a...)

Return the Kahan summation of a sequence of numbers `a::Vararg`. 

For more details about this method, see [Wikipedia](https://en.wikipedia.org/wiki/Kahan_summation_algorithm).

# Arguments
- `a::Vararg` Sequence of numbers to be summed.

# Returns
- The sum of the input numbers computed using Kahan summation. The return type is
  `eltype(a)'

# Error
Throw an `ArgumentError` when the argument `a` is empty.

# Example
```julia
KahanSum(1.0, 2.0, 3.0)
6.0

KahanSum(2//3, 6.7, BigFloat(5.6))
12.9666666666666667850904559600166976451873779296875
"""
function KahanSum(a::Vararg)
    !isempty(a) || throw(ArgumentError("The function KahanSum requires at least one argument"))
    T = eltype(a)
    sum::T = zero(T)
    c::T = zero(T)
    for x in a
        y::T = x - c
        t::T = sum + y
        c = (t - sum) - y
        sum = t
    end
    sum
end

# This function optionally uses a polylog2 function identity before
# it calls polylog2_helper. The polylog2 function has functional 
# relations for x --> 1/x,  x --> 1/(1-x), x --> (x-1)/x, and x --> 1-x. 
# But the linear convergence rate, given by cnd, is the same for x --> 1/x & 
# x --> 1/(1-x) and the same for x --> (x-1)/x & x --> 1-x. So we only choose
# between using x --> x, x --> 1/x, and x --> 1-x.

# The linear convergence rate is bounded above by 1/sqrt(3). The
# linear convergence rate for cis(pi/3) is 1/sqrt(3).

"""
    polylog2(x::Number)

Compute the numeric value of the dilogarithm. For a definition of this function,
see http://dlmf.nist.gov/25.12.E1. 

Examples:
```
julia> polylog2(0.125)
julia> 0.1291398601099534

The input can be a complex number:
julia>  polylog2(0.5 + 128.0*im)
julia> -12.176241112881845 + 7.648328923834815im

Finally, the input can be a BigFloat:

julia> setprecision(BigFloat,128);
julia> polylog2(BigFloat(0.125))
julia> 0.1291398601099534056689353043446094486239
```
"""
function polylog2(x::Number)
    T = typeof(x)
    cnd = x -> if isapprox(2, x, atol=eps(T)) Inf else abs2(x / (2 - x)) end
    c0 = cnd(x)
    c1 = if isapprox(0, x, atol=eps(T)) Inf else cnd(inv(x)) end
    c2 = cnd(1 - x)
    cmin = min(c0, c1, c2)

    R = if x == 0
        convert(T, 0), true
    elseif x == 1
        zeta2(T), true
    elseif cmin == c0 #no transformation
        q0 = x / (one(T) - x / 2)
        polylog2_helper(q0, x)
    elseif cmin == c1 #do x -> 1/x transformation
        q0 = inv(x - one(T)/2)
        f = polylog2_helper(q0, inv(x))
        #was -((f[1] + zeta2(T)) + clog(-x)^2 / 2), f[2]
        -KahanSum(zeta2(T), f[1], clog(-x)^2/2), f[2]
    else #do x -> 1-x transformation
        q0 = 2 * ((one(T) - x) / (one(T) + x))
        f = polylog2_helper(q0, one(T) - x)
        # I've experimented with replacing clog(one(T) - x))
        # with log1p(-x). It's not a clear win.
        #was  zeta2(T) - (f[1] + clog(x)*log1p(one(T)-x)), f[2]
        KahanSum(zeta2(T), -f[1], -clog(x)*clog(one(T)-x)), f[2]        
    end
    if R[2]
        R[1]
    else
        # When the running error bound is too great, we should
        # try again with a BigFloat with greater precision.
        error("Unable to evaluate(polylog2(", x, ")")
    end
end

"""
    mapabs(x)
When `x` is real, return abs(x); and when `x` is complex, return
abs(real(x)) + abs(real(x))im.
"""
function mapabs(x::Complex) 
    abs(real(x)) + im*abs(imag(x))
end

function mapabs(x::Real) 
    abs(x)
end

# return value of polylog(2,x) and a boolean that indicates success or failure.

# We have h = L + c (-x/2)^k. We could exploit this fact to extrapolate the limit
# and return early.

# It's a fun game to attempt to find the fewest number of Int64 add and multiplies 
# to compute (k+1)(k+2), (k+2)^2, and (k+3)*(k+4). Let's just let it be.

# We could use muladd to evaluate the dot product p0q0+p1q1+p2q2, but I'm not sure 
# we win. Julia's fma function doesn't allow complex arguments, so I'm not sure we
# gain any accuracy by using muladd.

# The value of he is a running error bound on the rounding error of h. For a 
# description of the running error see _Accuracy and Stability of Numerical Algorithms_,
# by Nicholas Higham (SIAM, 2002, ISBN 0-89871-521-0). 

function polylog2_helper(q0::Number, x::Number)
    T = typeof(x)
    #was q0 = x/(1-x/2)
    q1 = (-q0^2) / 4 # was -x^2/(4*(1-x/2)^2)
    q2 = (q0^3) / 9  # was x^3/(9*(1-x/2)^3)
    #was h = q0 + (q1 + q2) # not sure of best order to sum.
    h = KahanSum(q2,q1,q0) #q0 - q0^2/4 + q0^3/9
    N = 2^24 # magic number--it is a power of two for no particular reason
    k = zero(N)
    streak = zero(N)
    ks = zero(T) #Kahan summation corrector
    s0 = x/(x - 2)
    s1 = s0^2
    s2 = s0^3
    ε = eps(T)
    he = zero(T) # upper limit for error in h
    while k < N && streak < 5 && !isnan(h) && !isinf(h)  #magic number 5
        #was q3 = (-(k+1)*(k+2)*q0*x^3+(k+2)^2*q1*(x-2)*x^2+(k+3)*(k+4)*q2*(x-2)^2*x)/((k+4)^2*(x-2)^3)

        # We need to be careful with contagion. But these do Int64*float, and I think
        # these do the proper contagion.

        # It's fun to attempt to compute p0, p1, p2, and q3 with the least amount of
        # arithmetic. But the speedup is small, and the impact on accuracy is questionable.
        p0 = -((k+1)*(k+2))*s2
        p1 = (k+2)^2*s1
        p2 = ((k+3)*(k+4))*s0
        #was q3 = (p0 * q0 + (p1 * q1 + p2 * q2))/(k+4)^2 # not sure of best order to sum.
        q3 = KahanSum(p1*q1, p2*q2, p0*q0)/(k+4)^2
        qq3 = q3 - ks #start Kahan summation
        t = h + qq3
        ks = (t - h) - qq3
        streak = if (h == t) || isapprox(h,t,atol=ε) streak + 1 else 0 end
        h = t #end Kahan summation	
        he +=  mapabs(h)
        (q0,q1,q2) = (q1,q2,q3)
        k += 1
    end
    h, k < N && !isnan(h) && !isinf(h) && real(he) < 256*(1 + abs(real(h))) && imag(he) < 256*(1 + abs(imag(h)))
end

function polylog2(x::Int64)
    polylog2(convert(Float64,x)) 
end

function polylog2(x::Complex{Int64})
    polylog2(convert(Complex{Float64},x)) 
end

# polylog2(pi) = polylog2(convert(Float64,pi))
function polylog2(x::Irrational)     
    polylog2(convert(Float64,x))
end

# polylog2(im) = polylog2(convert(Complex{Float64},im))
function polylog2(x::Complex{Bool})     
    polylog2(convert(Complex{Float64},x))
end

function convergence_rate(x::Number)
    if isreal(x)
        α = -x/2
        μ = α/(1+α) # linear convergence rate
    else
      s = sqrt((x-1)/(conj(x)-1))
      α1 = x/(s-1)
      α2 = -x/(s+1)
      α = if abs2(α1/(α1 +1)) < abs2(α2/(α2+1)) α1 else α2 end
      μ = α/(α+1) # linear convergence rate
    end
    if isnan(μ) Inf else abs2(μ) end
end

function polylog2X(x::Number)
    T = typeof(x)
    μ0 = convergence_rate(x) # no transformation
    μ1 = convergence_rate(1/x) # x -> 1/x transformation
    μ2 = convergence_rate(1-x) # x -> 1-x transformation
    μ3 = convergence_rate(x/(x-1)) # x -> x/(x-1) transformation
    μmin = min(μ0, μ1, μ2, μ3)
    R = if x == 0
        convert(T, 0), true
    elseif x == 1
        zeta2(T), true
    elseif μmin == μ0 # no transformation
        #println("x -> x")
        polylog2X_helper(x)
    elseif μmin == μ1 # do x -> 1/x transformation
        # println("x -> 1/x")
        f = polylog2X_helper(1/x)
        #-((f[1] + zeta2(T)) + clog(-x)^2 / 2), f[2]
        -KahanSum(f[1], zeta2(T), clog(-x)^2 / 2), f[2]
    elseif μmin == μ2 # do x -> 1-x transformation 
        #println("x -> 1-x")
        f = polylog2X_helper(one(T) - x)
        # I've experimented with replacing clog(one(T) - x))
        # with log1p(-x). It's not a clear win.
        #zeta2(T) - (f[1] + clog(x)*clog(one(T)-x)), f[2]
        KahanSum(zeta2(T), -f[1], -clog(x)*clog(one(T)-x)), f[2]
    else # do x -> x -> x/(x-1) transformation
        #println("x -> x/(1-x)")
        f = polylog2X_helper(x/(x-1))
        -f[1]  - clog(1-x)^2/2, f[2]  # http://dlmf.nist.gov/25.12.E3
    end
    if R[2]
        R[1]
    else
        # When the running error bound is too large, we should
        # try again with a BigFloat with greater precision.
        error("Unable to evaluate(polylog2(", x, ")")
    end
end

# This code is based on a method that has a better linear convergence rate than 
# does polylog2_helper. Most experiments show that this code is _slower_ than
# polylog2 and no more accurate.
function polylog2X_helper(x)
    T = typeof(x)  
    if isreal(x)
        x = real(x)
        α = -x/2
        μ = α/(1+α) # linear convergence rate
    else
      s = sqrt((x-1)/(conj(x)-1))
      α1 = x/(s-1)
      α2 = -x/(s+1)
      α = if abs2(α1/(α1 +1)) < abs2(α2/(α2+1)) α1 else α2 end # not sure!
      μ = α/(α+1) # linear convergence rate
    end
    #println("|μ| = ", abs(μ))
    ks = zero(T) #Kahan summation corrector
    ε = eps(T)
    q0 = x/(1+α)
    q1 = x*(α+x/4)/(one(T)+α*(α+2)) #was: (x*α+x^2/4)/(α+1)^2 
    q2 = x*(α^2 + x*(α/2 + x/9))/((α+1)^3) #was: (x*α^2+(x^2*α)/2+x^3/9)/(α+1)^3
    N = 2^12
    k = zero(N)
    streak = zero(N)
    ks = zero(T) #Kahan summation corrector
    h = KahanSum(q0, q1, q2)
    he = zero(T) # running error

    #hoist some constants
    K1 = α^2*(α+x)/(α+1)^3

    K2 = 3*α+2*x
    K3 = 8*α+5*x
    K4 = α/(α+1)^2
    K2 *= K4
    K3 *= K4

    K5 = 3*α + x
    K6 = 10*α+3*x
    K7 = one(T)+α
    K5 /= K7
    K6 /= K7

    while k < N && streak < 5 && !isnan(h) && !isinf(h)
      #p0 = ((k+1)*(k+2)*α^2*(α+x))/((k+4)^2*(α+1)^3)
      #p1 = -((k+2)*α*(3*k*α+8*α+2*k*x+5*x))/((k+4)^2*(α+1)^2)
      #p2 = ((k+3)*(3*k*α+10*α+k*x+3*x))/((k+4)^2*(α+1))

      p0 = ((k+1)*(k+2))*K1
      p1 = -(k+2)*(K2*k + K3)
      p2 = (k+3)*(K5*k + K6)
      
      q3 = KahanSum(p1*q1, p2*q2, p0*q0)/(k+4)^2
      qq3 = q3 - ks #start Kahan summation
      t = h + qq3
      ks = (t - h) - qq3
      streak = if (h == t) || isapprox(h,t,atol=ε) streak + 1 else 0 end
      h = t #end Kahan summation	
      he +=  mapabs(h) #update running error
      (q0,q1,q2) = (q1,q2,q3)
      k += 1
    end
    #@show(k)
    #@show(he*eps(T))
    #@show(h)
    h, k < N && !isnan(h) && !isinf(h) && real(he) < 256*(1 + abs(real(h))) && imag(he) < 256*(1 + abs(imag(h)))
end
