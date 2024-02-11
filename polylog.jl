# Barton Willis, 2019, 2024

# This work is licensed under the CC0 1.0 Universal license.

# Julia code for the evaluation of the dilogarithm. The method is
# based on "The binomial transform of p-recursive sequences and 
# dilogarithm function," by Stephanie Harshbarger and Barton Willis.
# https://arxiv.org/pdf/1910.06928.pdf

# Extend eps.
import Base.eps
import Base.BigFloat
eps(::Type{Complex{T}}) where T <: AbstractFloat = eps(T)

"""
  clog(x::Number)

When x in (0,infinity) return log(x); otherwise, convert x to 
a complex and then dispatch log on x. This function is _not_
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
	@assert 3*cmin <= 1
    R = if x == 0
        convert(T, 0), true
    elseif x == 1
        convert(T, pi)^2 / 6, true
    elseif cmin == c0 #no transformation
        q0 = x / (one(T) - x / 2)
        polylog2_helper(q0, x)
    elseif cmin == c1 #do x -> 1/x transformation
        q0 = inv(x - one(T)/2)
        f = polylog2_helper(q0, inv(x))
        -f[1] - convert(T, pi)^2 / 6 - clog(-x)^2 / 2, f[2]
    else #do x -> 1-x transformation
        q0 = 2 * ((one(T) - x) / (one(T) + x))
        f = polylog2_helper(q0, one(T) - x)
        # I don't think chaning clog(1-x) to log1p(-x) is a win?
        -f[1] + convert(T, pi)^2 / 6 - clog(x) * clog(one(T) - x), f[2]
    end
    if R[2]
        R[1]
    else
        error("Unable to evaluate(polylog2(", x, ")")
    end
end

# return value of polylog(2,x) and a boolean that indicates success or failure.

# We have h = L + c (-x/2)^k. We could exploit this fact to extrapolate the limit
# and return early.

# It's a fun game to attempt to find the fewest number of Int64 add and multiplies 
# to compute (k+1)(k+2), (k+2)^2, and (k+3)*(k+4). Let's just let it be.

# We could use muladd to evaluate the dotproduct p0q0+p1q1+p2q2, but I'm not sure 
# we win. Julia's fma function doesn't allow complex arguments, so I'm not sure we
# gain any accuracy by using muladd

function polylog2_helper(q0::Number, x::Number)
    T = typeof(x)
    #was q0 = x/(1-x/2)
    q1 = (-q0^2) / 4 # was -x^2/(4*(1-x/2)^2)
    q2 = (q0^3) / 9  # was x^3/(9*(1-x/2)^3)
    h = q0 + (q1 + q2) # not sure of best order to sum.
    N = 2^24 # magic number--it is a power of two for no particular reason
    k = zero(N)
    streak = zero(N)
    cndR = abs(real(q0)) + abs(real(q1)) + abs(real(q2)) #real part of sum condition number.
    cndI = abs(imag(q0)) + abs(imag(q1)) + abs(imag(q2)) #imaginary part of sum condition number.
    ks = zero(T) #Kahan summation corrector
    s0 = x / (x - 2)
    s1 = s0^2
    s2 = s0^3
    ε = eps(T)
    while k < N && streak < 5 && !isnan(h) && !isinf(h)  #magic number 5
        #was q3 = (-(k+1)*(k+2)*q0*x^3+(k+2)^2*q1*(x-2)*x^2+(k+3)*(k+4)*q2*(x-2)^2*x)/((k+4)^2*(x-2)^3)

        # We need to be careful with contagion. But these do Int64*float, and I think
        # these do the proper contagion.
        p0 = -(k+1)*(k+2)*s2
        p1 = (k+2)^2*s1
        p2 = (k+3)*(k+4)*s0
        q3 = (p0 * q0 + (p1 * q1 + p2 * q2))/(k+4)^2 # not sure of best order to sum.
        qq3 = q3 - ks #start Kahan summation
        t = h + qq3
        ks = (t - h) - qq3
        streak = if (h == t) || isapprox(h,t,atol=ε) streak + 1 else 0 end
        h = t #end Kahan summation	
        cndR += abs(real(q3))
        cndI += abs(imag(q3))
        (q0,q1,q2) = (q1,q2,q3)
        k += 1
    end
    h, k < N && !isnan(h) && !isinf(h) && cndR < 16 && cndI < 16
end

function polylog2(x::Int64)
    polylog2(convert(Float64,x)) 
end