# Test the identity http://dlmf.nist.gov/25.12.E3
function dlmf_25_12_3(x)
   isapprox(polylog2(x) + polylog2(x/(x-1)), -log(1-x)^2 / 2, rtol = 16*eps(x))
end

