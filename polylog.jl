# Barton Willis, Copyright 2019, 2024

# This work is licensed under the CC0 1.0 Universal license.

# Julia code for the evaluation of the dilogarithm. The method is
# based on "The binomial transform of p-recursive sequences and 
# dilogarithm function," by Stephanie Harshbarger and Barton Willis.
# https://arxiv.org/pdf/1910.06928.pdf

# Extend eps.
import Base.eps
eps(x::Type{Complex{Float16}}) = eps(Float16)
eps(x::Type{Complex{Float32}}) = eps(Float32)
eps(x::Type{Complex{Float64}}) = eps(Float64)
eps(x::Type{Complex{BigFloat}}) = eps(BigFloat)

# convert a complex to a bigfloat
BigFloat(x::Complex) = BigFloat(real(x)) + BigFloat(imag(x))im

"""
  mylog(x::Number)

When x in (0,infinity) return log(x); otherwise, convert x to 
a complex and then dispatch log on x. This function is _not_
intended to be a user-level function.
"""
function mylog(x::Number)
	if isreal(x) && x > 0
		log(x)
	else 
		log(Complex(x))	
	end
end

# This function optionally uses a polylog2 function identity before
# it calls polylog2_helper. The polylog2 function has functional 
# relations for x --> 1/x,  x --> 1/(1-x), x --> (x-1)/x, and x --> 1-x. 
# But the convergence rate, given by cnd,  is the same for x --> 1/x & 
# x --> 1/(1-x) and the same for x --> (x-1)/x & x --> 1-x. So we only choose
# betwen using x --> x, x --> 1/x, and x --> 1-x.

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
	cnd = x -> if isapprox(2,x,atol=eps(T)) Inf else abs2(x/(2-x)) end
	c0 = cnd(x)
	c1 = if isapprox(0,x,atol=eps(T)) Inf else cnd(1/x) end 
	c2 = cnd(1-x)
	cmin = min(c0,c1,c2)
	R = if x == 0
		   convert(T,0), true
	    elseif x == 1
		   convert(T, pi)^2/6, true
	    elseif cmin == c0 #no transformation
		   q0 = x/(1-x/2)
		   polylog2_helper(q0,x)
	    elseif cmin == c1 #do x -> 1/x transformation
		   q0 = 1/(x-1//2)
		   f = polylog2_helper(q0,1/x)
		  -f[1] - convert(T, pi)^2/6 - mylog(-x)^2/2,f[2]
	    else #do x -> 1-x transformation
		   q0 = 2*((1-x)/(1+x))
		   f = polylog2_helper(q0,1-x)
		  -f[1] + convert(T, pi)^2/6 - mylog(x)*mylog(1-x),f[2]
	    end
	if R[2]
		R[1]
	else 
		error("Unable to evaluate(polylog2(", x, ")")
	end
end

# return value of polylog(2,x), the condition number of the sum, the 
# number of terms summed, and a boolean that indicates sucess or failure.
function polylog2_helper(q0::Number, x::Number)
	T = typeof(x)
	#was q0 = x/(1-x/2)
    q1 = -q0^2/4 # was -x^2/(4*(1-x/2)^2)
	q2 = q0^3/9  # was x^3/(9*(1-x/2)^3)
	h = q0+(q1+q2)
	N = convert(Int64, 2^24) #magic number--it is a power of two for no particular reason
    k = zero(N)
    streak = zero(N)
	cndR = abs(real(q0))+abs(real(q1))+abs(real(q2)) #real part of sum condition number.
	cndI = abs(imag(q0))+abs(imag(q1))+abs(imag(q2)) #imaginary part of sum condition number.
	ep = eps(T)
	ks = zero(T) #Kahan summation corrector
	s0 = x/(x-2)
	s1 = s0^2
	s2 = s0^3

    while k < N && streak < 5 && !isnan(h) && !isinf(h)  #magic number 5
      #was q3 = (-(k+1)*(k+2)*q0*x^3+(k+2)^2*q1*(x-2)*x^2+(k+3)*(k+4)*q2*(x-2)^2*x)/((k+4)^2*(x-2)^3)
	  
	  # We need to be careful with contagion. Replacing ((k+3)*s0)/(k+4)
	  # by ((k+3)/(k+4))*s0 is OK when s0 is a binary64, but not OK when s0 is
	  # a BigFloat. So we do (integer x Float)/integer, and I think this is OK.
	  p0 = -((k+1)*(k+2)*s2)/((k+4)^2)
	  p1 = ((k+2)^2*s1)/((k+4)^2)
	  p2 = ((k+3)*s0)/(k+4)
	  q3 = p0*q0 + (p1*q1 + p2*q2) # not sure of best order to sum.
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
	OK = k < N && !isnan(h) && !isinf(h) && cndR < 16 && cndI < 16
	h, OK
end