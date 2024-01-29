# Barton Willis, Copyright 2019, 2024

# This work is licensed under the CC0 1.0 Universal license.

# Julia code for the evaluation of the dilogarithm. The method is
# based on "The binomial transform of p-recursive sequences and 
# dilogarithm function," by Stephanie Harshbarger and Barton Willis.
# https://arxiv.org/pdf/1910.06928.pdf

# Extend eps.
import Base.eps
eps(x::Int) = zero(x)
eps(x::Rational) = zero(x)
eps(x::Rational{Int64}) = zero(x)
eps(x::Type{Rational{Int64}}) = zero(x)
eps(x::BigInt) = zero(x)
eps(x::Complex{Int64}) = zero(Int64)
eps(x::Complex{Rational{Int64}}) = zero(Int64)
eps(::Complex{Rational{BigInt}}) = zero(Rational{BigInt})
eps(x::Complex{Float16}) = eps(Float16)
eps(x::Complex{Float32}) = eps(Float32)
eps(x::Complex{Float64}) = eps(Float64)
eps(x::Type{Complex{Float64}}) = eps(Float64)
eps(x::Complex{BigFloat}) = eps(BigFloat)
eps(x::Type{Complex{BigFloat}}) = eps(BigFloat)

# convert a complex to a bigfloat
BigFloat(x::Complex) = BigFloat(real(x)) + BigFloat(imag(x))im

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
	# call polylog2_transform & check for success
	f = polylog2_transform(x::Number)
	if f[4]
		f[1]
	else
		error("Unable to evaluate(polylog2(", x, ")")
	end 
end

# When x in (0,infinity) return log(x); otherwise, convert x to 
# a complex and then dispatch log. This helps prevent some X + 0.0im
# results.
function mylog(x::Number)
	if isreal(x) && x > 0
		log(x)
	else 
		log(Complex(x))	
	end
end

# This function optionally uses a polylog2 function identity before
# it calls polylog2_helper. This is not intended to be a user level 
# function. Presumably, it chooses the identity to gain speed and
# accuracy.

# The polylog2 function has functional relations for x --> 1/x, 
# x --> 1/(1-x), x --> (x-1)/x, and x --> 1-x. But the convergence
# rate, given by cnd,  is the same for x --> 1/x & x --> 1/(1-x)
# and the same for x --> (x-1)/x & x --> 1-x. So we only choose
# betwen using x --> x, x --> 1/x, and x --> 1-x.
function polylog2_transform(x::Number)
	T = typeof(x)
	cnd = x -> if x == 2 Inf else abs2(x/(2-x)) end
	c0 = cnd(x)
	c1 = if x == 0 Inf else cnd(1/x) end 
	c2 = cnd(1-x)
	cmin = min(c0,c1,c2)
	if x == 0
		convert(T,0), convert(T,0), 0, true
	elseif x == 1
		convert(T, pi^2/6), convert(T,0), 0, true
	elseif cmin == c0 #no transformation
		q0 = x/(1-x/2)
		polylog2_helper(q0,x)
	elseif cmin == c1 #do x -> 1/x transformation
		q0 = 1/(x-1/2)
		f = polylog2_helper(q0,1/x)
		-f[1] - convert(T, pi^2/6) - mylog(-x)^2/2,f[2], f[3], f[4]
	else #do x -> 1-x transformation
		q0 = 2*((1-x)/(1+x))
		f = polylog2_helper(q0,1-x)
		-f[1] + convert(T, pi^2/6) - mylog(x)*mylog(1-x),f[2], f[3], f[4]
	end
end
	
# return value of polylog(2,x), the condition number of the sum, the 
# number of terms summed, and a boolean that indicates sucess or failure.
function polylog2_helper(q0::Number, x::Number)
	T = typeof(x)
	#was q0 = x/(1-x/2)
    q1 = -q0^2/4 # was -x^2/(4*(1-x/2)^2)
	q2 = q0^3/9  # was x^3/(9*(1-x/2)^3)
	h = @evalpoly(q0, 0, 1,-1/4, 1/9) #was h = q0+(q1+q2)
	N = convert(Int64, 2^24) #magic number--it is a power of two for no particular reason
    k = zero(N)
    streak = zero(N)
	cndR = abs(real(q0))+abs(real(q1))+abs(real(q2)) #real part of sum condition number.
	cndI = abs(imag(q0))+abs(imag(q1))+abs(imag(q2)) #imaginary part of sum condition number.
	ep = eps(T)
	ks = zero(T) #Kahan summation corrector
    w0 = x/(x-2)^3 #hoist two constants by hand
	w1 = (x-2)^2
    while k < N && streak < 5 #magic number 5
      #was q3 = (-(k+1)*(k+2)*q0*x^3+(k+2)^2*q1*(x-2)*x^2+(k+3)*(k+4)*q2*(x-2)^2*x)/((k+4)^2*(x-2)^3)
	  q3 = w0*((x*(k+2)* ((k+2)*q1*(x-2)-(k+1)*q0*x)) +w1*(k+3)*(k+4)*q2)/((k+4)^2)
	  qq3 = q3-ks #start Kahan summation
	  t = h+qq3 
	  ks = (t - h) - qq3
	  streak = if (h == t) streak + 1 else 0 end
	  h = t #end Kahan summation	
	  cndR += abs(real(q3))
	  cndI += abs(imag(q3))	
	  q0 = q1
      q1 = q2
	  q2 = q3	
      k += 1
	end
	h, ep*(if cndI == 0 cndR else cndR + cndI*im end), k, k < N && !isnan(h) && !isinf(h)
end
	

