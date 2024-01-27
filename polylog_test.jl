# Test the identity http://dlmf.nist.gov/25.12.E3 .
# This identity is valid off [1,infty). For a input in [1,infy),
# return true.
function dlmf_25_12_3(x)
   if imag(x)==0 && x >= 1
      true
   else
     isapprox(polylog2(x) + polylog2(x/(x-1)), -log(1-x)^2 / 2, rtol = 16*eps(x))
end

