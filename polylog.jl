# Barton Willis, Copyright 2019, 2024

# This work is licensed under the CC0 1.0 Universal license
# Julia code for the evaluation of the dilogarithm. The method is
# based on "The binomial transform of p-recursive sequences and 
# dilogarithm function," Stephanie Harshbarger and Barton Willis.
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
    polylog2(x)

Compute the numeric value of the dilogarithm. For a definition of this function,
see http://dlmf.nist.gov/25.12.E1 . The function returns a four tuple of the 
for (dilogarithm(x), cnd, n, bool), where cnd is the condition sum that is used
for evalution, n is the number of terms summed, and bool is a boolean that 
indicates sucess or failure.

```
julia> polylog2(1.01)
julia> (1.700732144324037 - 0.03125988630910074im, 2.2149194678541363e-18, 9, true)

The input can be a complex number:
julia> polylog2(3.1 + 5.1*im)
julia> -1.0831276752057937 + 3.9312546561253905im, 1.8745409421641466e-17 + 3.595690852077723e-17im, 17, true)

Finally, the input can be a BigFloat:

julia> polylog2(BigFloat(1.01))
julia> (1.700732144324037086378889112827585843561092568405653427462945581660835230870607 - 0.03125988630910073633439461464037805326512443369418564280068429886030330779614308im, 1.722934710961044873101367377667284112099228890489560282524331651371580237533837e-79, 36, true)

```
"""
function polylog2(x::Number)
	T = typeof(x)
	cnd = x -> if x == 2 Inf else abs2(x/(2-x)) end
	c0 = cnd(x)
	c1 = if abs2(x) <= 1 Inf else cnd(1/x) end #only do 1/x transform for x outside unit cirle.
	c2 = cnd(1-x)
	cmin = min(c0,c1,c2)
	if x == 0
		convert(T,0), convert(T,0), 0, true
	elseif x == 1
		convert(T, pi^2/6), convert(T,0), 0, true
	elseif cmin == c0 #no transformation
		polylog2_helper(x)
	elseif cmin == c1 #do x -> 1/x transformation
		f = polylog2_helper(1/x)
		-f[1] - convert(T, pi^2/6) - log(Complex(-x))^2/2,f[2], f[3], f[4]
	else #do x -> 1-x transformation
		f = polylog2_helper(1-x)
		-f[1] + convert(T, pi^2/6) - log(Complex(x))*log(Complex(1-x)),f[2], f[3], f[4]
	end
end
	
# return value of polylog(2,x), the condition number of the sum, the 
# number of terms summed, and a boolean that indicates sucess or failure.
function polylog2_helper(x::Number)
	T = typeof(x)
	q0 = x/(1-x/2)
    q1 = -x^2/(4*(1-x/2)^2)
	q2 = x^3/(9*(1-x/2)^3)
    h = q0+(q1+q2)
	
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
      #q3 = (-(k+1)*(k+2)*q0*x^3+(k+2)^2*q1*(x-2)*x^2+(k+3)*(k+4)*q2*(x-2)^2*x)/((k+4)^2*(x-2)^3)
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
	

