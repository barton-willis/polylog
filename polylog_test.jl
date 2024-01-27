# Test the identity http://dlmf.nist.gov/25.12.E3 .
# This identity is valid off [1,infty). For a input in [1,infy),
# return true.
function dlmf_25_12_3_E3(x)
   if imag(x)==0 && x >= 1
      true
   else
     isapprox(polylog2(x) + polylog2(x/(x-1)), -log(1-x)^2 / 2, rtol = 16*eps(x))
end

# http://dlmf.nist.gov/25.12.E5
function dlmf_25_12_E5(z,m::Int64)
   if abs(z) < 1 && m > 0 
      s = 0
      k = 0
      while k < m
         s += polylog2(z*cis(2*pi*k/m))
         k += 1
      end
      isapprox(polylog2(z^m), m*s, rtol = 16*eps(z))
   else
      true
   end
end
