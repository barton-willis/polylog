# Many of these tests use function identities that the code uses to 
# transform the input to gain better convergence. Such tests are not
# particularly good accuracy tests.

# A modified relative difference function. When both |a| >= 1 and |b| >= 1 and
# a & b are real, return abs(a-b)/(ε * min(|a|, |b|))

#Maybe this needs to be modified for denormal numbers?
function rd(a::Real,b::Real)
   if isinf(a) || isinf(b) || isnan(a) || isnan(b)
      Inf
   else
    floor(Int64, abs(a-b)/(eps(typeof(a)) * max(1, min(abs(a),abs(b)))))
   end
end

function rd(a::Number,b::Number)
     max(rd(real(a),real(b)), rd(imag(a),imag(b)))
end

"""
    boa(a::Real, b::Real)

Calculate the number of bits of agreement between two real numbers `a` and `b`.
Specifically, when `a` and `b` are both nonzero and have the same sign, return 

     floor(Int64, min(precision(a), -log2(|a-b|/min(|a|,|b|))

When `a` and `b` have opposite signs, return zero; and when either `a` or `b` is
zero, return  floor(Int64, min(precsion(a), -log2(max(|a|,|b|)))).

Also, when `a` and `b` don't have the same type, we apply `promote` to them.

# Examples
```julia
julia> boa(2.0, 2.0+ 1/2^10)
11

julia> boa(0.0, 0.0)
53

julia> boa(-5.0, 5.0)
0
"""
function boa(a::Real, b::Real)
   (a,b) = promote(a,b)
   m =  precision(a)
   q =  if iszero(a) && iszero(b) # both a & b are zero
           m
         elseif sign(a) != sign(b) 
           0
         elseif iszero(a)
           -log2(abs(b))
         elseif iszero(b)
            -log2(abs(a))
         else
            -log2(abs((a-b)/min(abs(a),abs(b))))
         end
   max(0, floor(Int64, min(m,q)))
 end

function boa(x::Complex,y::Complex)
   min(boa(real(x),real(y)),boa(imag(x),imag(y)))
end

function test_report(results)
   Q = sort(collect(results), by=x->x[1], rev=true)
   s = 0
   n = 0
   w = typemax(Int64)
   for x in Q 
      n += x[2]
      s += x[1]*x[2]
      w = min(w,x[1])
   end
   println("number of tests = ", n)
   println("average correct bits = ",s/n)
   println("worst correct bits = ",w)
   Q
end

# Test the identity http://dlmf.nist.gov/25.12.E3 . This identity is valid 
# off [1,infinity), but this function does not check this condition.
function dlmf_25_12_3_E3(x)
   #rd(polylog2(x) + polylog2(x/(x-1)), -log(1-x)^2 / 2)
   boa(polylog2(x) + polylog2(x/(x-1)), -log(1-x)^2 / 2)
end

# Test dlmf_25_12_3_E3 inside and on the unit circle.
function polylog2_test1(T::DataType,m::Int64,n::Int64)
   results = Dict()
   for i = 0 : m
      for j = 0 : n
         x = convert(T,(i/m) * cis(2*pi*j/n))
         if !iszero(imag(x)) || real(x) < 1
            q = dlmf_25_12_3_E3(x)
            results[q] = if haskey(results,q) results[q]+1 else 1 end
         end
      end
   end
   test_report(results)
end

# Test dlmf_25_12_3_E3 on the unit circle
function polylog2_test2(T::DataType,n::Int64)
   results = Dict()
   for j = 1 : n-1
       x = convert(T, cis(2*pi*j/n))
       q = dlmf_25_12_3_E3(x)
       results[q] = if haskey(results,q) results[q]+1 else 1 end
   end
   test_report(results)
end

# http://dlmf.nist.gov/25.12.E5
function dlmf_25_12_E5(T,z,m::Int64)
   if abs(z) < 1 && m > 0
      s = [zero(T)]
      k = 0
      while k < m
         θ = (2*k*convert(T,pi))/m
         push!(s, polylog2(z*cis(θ)))
         k += 1
      end
      s = KahanSum(s...)
      boa(polylog2(z^m), m*s)
   else
      true
   end
end

function polylog2_test3(T::DataType,n::Int64,m::Int64)
   results = Dict()
   for i = 1 : n-1
      for j = 0 : n
         x = convert(T, (i/n)*cis(2*pi*j/n))
         q = dlmf_25_12_E5(T,x,m)
         results[q] = if haskey(results,q) results[q]+1 else 1 end
      end
   end
   test_report(results)
end

#http://dlmf.nist.gov/25.12.E7
function dlmf_25_12_E7(T,x)
   if zero(typeof(x)) <= x && x <= convert(typeof(x), 2*pi)
     pie = convert(T,pi)
     boa(real(polylog2(cis(x))),  pie*(pie/6 - x/2) + x^2/4)
   else 
      true
   end
end

function polylog2_test4(T::DataType, n::Int64)
   results = Dict()
   pie = convert(T,pi)
   for i = 0 : n-1
       q = dlmf_25_12_E7(T,(2*pie*i)/n)
       results[q] = if haskey(results,q) results[q]+1 else 1 end
   end
   test_report(results)
end

#See https://en.wikipedia.org/wiki/Dilogarithm
function polylog2_id_1(x)
   boa(polylog2(x)+polylog2(-x),polylog2(x^2)/2)
end

function polylog2_test5(T::DataType,n::Int64)
   results = Dict()
   pie = convert(T, pi)
   for i = 1 : n
      for j = 0 : n-1
         x = ((convert(T,2)*i)/n) * cis((2*pie*j)/n)
         q =  polylog2_id_1(x)
         results[q] = if haskey(results,q) results[q]+1 else 1 end
      end
   end
   test_report(results)
end