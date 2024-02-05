# Many of these tests use function identities that the code uses to 
# transform the input to gain better convergence. Such tests are not particularly
# good.

# A modified relative difference function. Maybe this needs to be modified
# for denormal numbers?
function rd(a,b)
    floor(Int64, abs(a-b)/(eps(typeof(a)) * max(1, min(abs(a),abs(b)))))
end

function test_report(results)
   Q = sort(collect(results), by=x->x[1])
   s = 0
   n = 0
   w = 0
   for x in Q 
      n += x[2]
      s += x[1]*x[2]
      w = max(w,x[1])
   end
   println("number of tests = ", n)
   println("average error = ",s/n)
   println("worst error = ",w)
   Q
end

# Test the identity http://dlmf.nist.gov/25.12.E3 .
# This identity is valid off [1,infinity). For a input in [1,infinity),
# return true.
function dlmf_25_12_3_E3(x)
   if imag(x)==0 && real(x) >= 1
      true
   else
     rd(polylog2(x) + polylog2(x/(x-1)), -log(1-x)^2 / 2)
   end
end

# Test dlmf_25_12_3_E3 inside the unit circle.
function polylog2_test1(T::DataType,n::Int64)
   results = Dict()
   while n > 0
       x = convert(T,rand() * cis(2*pi*rand()))
       q = dlmf_25_12_3_E3(x)
       results[q] = if haskey(results,q) results[q]+1 else 1 end
       n -= 1
   end
   test_report(results)
end

# Test dlmf_25_12_3_E3 on the unit circle
function polylog2_test2(T::DataType,n::Int64)
   results = Dict()
   while n > 0
       x = convert(T, cis(2*pi*rand()))
       q = dlmf_25_12_3_E3(x)
       results[q] = if haskey(results,q) results[q]+1 else 1 end
       n -= 1
   end
   test_report(results)
end

# http://dlmf.nist.gov/25.12.E5
function dlmf_25_12_E5(T,z,m::Int64)
   if abs(z) < 1 && m > 0 
      s = zero(T)
      k = 0
      while k < m
         θ = (2*k*convert(T,pi))/m
         s += polylog2(z*cis(θ))
         k += 1
      end
      rd(polylog2(z^m), m*s)
   else
      true
   end
end

function polylog2_test3(T::DataType,n::Int64)
   results = Dict()
   while n > 0
       x = convert(T, rand()*cis(2*pi*rand()))
       q = dlmf_25_12_E5(T,x,rand(1:4))
       results[q] = if haskey(results,q) results[q]+1 else 1 end
       n -= 1
   end
   test_report(results)
end

#http://dlmf.nist.gov/25.12.E7
function dlmf_25_12_E7(T,x)
   if zero(typeof(x)) <= x && x <= convert(typeof(x), 2*pi)
     pie = convert(T,pi)
     rd(real(polylog2(cis(x))),  pie*(pie/6 - x/2) + x^2/4)
   else 
      true
   end
end

function polylog2_test4(T::DataType, n::Int64)
   results = Dict()
   while n > 0
       x = convert(T,2*pi*rand())
       q =  dlmf_25_12_E7(T,x)
       results[q] = if haskey(results,q) results[q]+1 else 1 end
       n -= 1
   end
   test_report(results)
end

#See https://en.wikipedia.org/wiki/Dilogarithm
function polylog2_id_1(x)
   rd(polylog2(x)+polylog2(-x),polylog2(x^2)/2)
end

function polylog2_test5(T::DataType,n::Int64)
   results = Dict()
   while n > 0
       x = convert(T, 2.0 * rand() * cis(2*pi*rand()))
       q =  polylog2_id_1(x)
       results[q] = if haskey(results,q) results[q]+1 else 1 end
       n -= 1
   end
   test_report(results)
end