using Test

function myprintln(s::String)
    printstyled(s, bold=true, color = :red)
    println()
end
φ = Base.MathConstants.golden

# Special values Int64 See https://en.wikipedia.org/wiki/Dilogarithm

println()
myprintln("Special Values Test")
@testset begin  
    @test polylog2(0) == 0
    @test polylog2(1) ≈ pi^2/6
    @test polylog2(2) ≈ pi^2/4 - im*pi*log(2)
end 

function polylog2_f16(x) 
    polylog2(convert(Float16,x))
end

pi16 = convert(Float16,pi)
φ16 = convert(Float16,φ)
ε = eps(Float16) 

function log16(x)
    log(convert(Float16,x))
end
# Special values using binary16 See https://en.wikipedia.org/wiki/Dilogarithm
println()
myprintln("Binary 16 Tests")
@testset begin  
    @test polylog2_f16(1/3)-polylog2_f16(1/9)/6 ≈ pi16^2/18 - log16(3)^2/6 atol = ε
    @test polylog2_f16(-1/3)-polylog2_f16(1/9)/3 ≈ -pi16^2/18 + log16(3)^2/6 atol = ε
    @test polylog2_f16(-1/2)+polylog2_f16(1/9)/6 ≈ -pi16^2/18 + log16(2)*log16(3)-log16(2)^2/2-log16(3)^2/3 atol = ε 
    @test polylog2_f16(1/4)+polylog2_f16(1/9)/3 ≈ pi16^2/18+2*log16(2)*log16(3)-2*log16(2)^2-(2/3)*log16(3)^2 atol = ε
    @test polylog2_f16(-1/8)+polylog2_f16(1/9) ≈ -log(9/8)^2/2 atol = ε 
    @test 36*polylog2_f16(1/2)-36*polylog2_f16(1/4)-12*polylog2_f16(1/8)+6*polylog2_f16(1/64) ≈ pi16^2 atol = ε
    @test polylog2_f16(-1.0) ≈ -pi16^2/12 atol = ε
    @test polylog2_f16(0.0) == 0.0 
    @test polylog2_f16(1/2) ≈ pi16^2/12 - log16(2)^2/2 atol = ε
    @test polylog2_f16(1.0) ≈ pi16^2/6 atol = ε
    @test polylog2_f16(-1/φ16) ≈ -pi16^2/15 + log(φ16)^2/2 atol = ε
    @test polylog2_f16(-φ16) ≈ -pi^2/10 - log(φ16)^2 atol = ε
    @test polylog2_f16(2-φ16) ≈ pi16^2/15 - log(φ16)^2 atol = ε 
    @test polylog2_f16(1/φ16) ≈ pi16^2/10 - log(φ16)^2 atol = ε 
    @test polylog2_f16(sqrt(2)-1)-polylog2_f16(1-sqrt(2)) ≈ pi16^2/8 - log16(1+sqrt(2))^2/2 atol = ε
    @test polylog2_f16(φ16) ≈ 11*pi16^2/15 + clog(-1/φ16)^2/2 atol = 8*ε
    @test polylog2_f16(φ16^2) ≈ -11*pi16^2/15 - clog(-φ16)^2 atol = 8*ε
end

function polylog2_f32(x) 
    polylog2(convert(Float32,x))
end

pi32 = convert(Float32,pi)
φ32 = convert(Float32,φ)

function log32(x)
    log(convert(Float32, x))
end
ε = eps(Float32)
# Special values using binary32 See https://en.wikipedia.org/wiki/Dilogarithm

println()
myprintln("Binary32 Tests")
@testset begin  
    @test polylog2_f32(1/3)-polylog2_f32(1/9)/6 ≈ pi32^2/18 - log32(3)^2/6 atol=ε
    @test polylog2_f32(-1/3)-polylog2_f32(1/9)/3 ≈ -pi32^2/18 + log32(3)^2/6 atol=ε
    @test polylog2_f32(-1/2)+polylog2_f32(1/9)/6 ≈ -pi32^2/18 + log32(2)*log32(3)-log32(2)^2/2-log32(3)^2/3 atol=ε
    @test polylog2_f32(1/4)+polylog2_f32(1/9)/3 ≈ pi32^2/18+2*log32(2)*log32(3)-2*log32(2)^2-(2/3)*log32(3)^2 atol=ε
    @test polylog2_f32(-1/8)+polylog2_f32(1/9) ≈ -log32(9/8)^2/2 atol=ε
    @test 36*polylog2_f32(1/2)-36*polylog2_f32(1/4)-12*polylog2_f32(1/8)+6*polylog2_f32(1/64) ≈ pi32^2 atol=ε
    @test polylog2_f32(-1.0) ≈ -pi32^2/12 atol=ε
    @test polylog2_f32(0.0) == 0.0f0 
    @test polylog2_f32(1/2) ≈ pi32^2/12 - log32(2)^2/2 atol=ε
    @test polylog2_f32(1.0) ≈ pi32^2/6 atol=ε
    @test polylog2_f32(-1/φ32) ≈ -pi32^2/15 + log32(φ32)^2/2 atol=ε
    @test polylog2_f32(-φ32) ≈ -pi32^2/10 - log32(φ32)^2 atol=ε
    @test polylog2_f32(2-φ32) ≈ pi32^2/15 - log32(φ32)^2 atol=ε
    @test polylog2_f32(1/φ32) ≈ pi32^2/10 - log32(φ32)^2 atol=ε
    @test polylog2_f32(sqrt(2)-1)-polylog2_f32(1-sqrt(2)) ≈ pi32^2/8 - log32(1+sqrt(2))^2/2 atol=ε
    @test polylog2_f32(φ32) ≈ 11*pi32^2/15 + clog(-1/φ32)^2/2 atol=ε
    @test polylog2_f32(φ32^2) ≈ -11*pi32^2/15 - clog(-φ32)^2 atol=4*ε
end

ε = eps(Float64)
# Special values using binary64 See https://en.wikipedia.org/wiki/Dilogarithm
println()
myprintln("Binary64 Tests")
@testset begin  
    @test polylog2(1/3)-polylog2(1/9)/6 ≈ pi^2/18 - log(3)^2/6 atol = ε
    @test polylog2(-1/3)-polylog2(1/9)/3 ≈ -pi^2/18 + log(3)^2/6 atol = ε
    @test polylog2(-1/2)+polylog2(1/9)/6 ≈ -pi^2/18 + log(2)*log(3)-log(2)^2/2-log(3)^2/3 atol = ε
    @test polylog2(1/4)+polylog2(1/9)/3 ≈ pi^2/18+2*log(2)*log(3)-2*log(2)^2-(2/3)*log(3)^2 atol = 4*ε
    @test polylog2(-1/8)+polylog2(1/9) ≈ -log(9/8)^2/2 atol = ε 
    @test 36*polylog2(1/2)-36*polylog2(1/4)-12*polylog2(1/8)+6*polylog2(1/64) ≈ pi^2 atol = 8*ε
    @test polylog2(-1.0) ≈ -pi^2/12 atol = ε
    @test polylog2(0.0) == 0.0   
    @test polylog2(1/2) ≈ pi^2/12 - log(2)^2/2 atol = ε
    @test polylog2(1.0) ≈ pi^2/6 atol = ε 
    @test polylog2(-1/φ) ≈ -pi^2/15 + log(φ)^2/2 atol = ε
    @test polylog2(-φ) ≈ -pi^2/10 - log(φ)^2 atol = ε
    @test polylog2(2-φ) ≈ pi^2/15 - log(φ)^2 atol = ε
    @test polylog2(1/φ) ≈ pi^2/10 - log(φ)^2 atol = ε 
    @test polylog2(sqrt(2)-1)-polylog2(1-sqrt(2)) ≈ pi^2/8 - log(1+sqrt(2))^2/2 atol = 4*ε
    @test polylog2(convert(Float64,φ)) ≈ 11*pi^2/15 + log(Complex(-1/φ))^2/2 atol = 4*ε
    @test polylog2(φ^2) ≈ -11*pi^2/15 - log(Complex(-φ))^2 atol = 4*ε
end

ε = eps(BigFloat)
piBF = convert(BigFloat,pi)
φBF = convert(BigFloat,φ)
function polylog2_bigfloat(x)
    polylog2(convert(BigFloat, x))
end

function logBF(x)
    clog(convert(BigFloat,x))
end

# Special values using binary64 See https://en.wikipedia.org/wiki/Dilogarithm
println()
myprintln("BigFloat Tests")
@testset begin  
    @test polylog2_bigfloat(1//3)-polylog2_bigfloat(1//9)/6 ≈ piBF^2/18 - logBF(3)^2/6 atol = ε
    @test polylog2_bigfloat(-1//3)-polylog2_bigfloat(1//9)/3 ≈ -piBF^2/18 + logBF(3)^2/6 atol = ε
    @test polylog2_bigfloat(-1//2)+polylog2_bigfloat(1//9)/6 ≈ -piBF^2/18 + logBF(2)*logBF(3)-logBF(2)^2/2-logBF(3)^2/3 atol = ε
    @test polylog2_bigfloat(1//4)+polylog2_bigfloat(1//9)/3 ≈ piBF^2/18+2*logBF(2)*logBF(3)-2*logBF(2)^2-(2//3)*logBF(3)^2 atol = 4*ε
    @test polylog2_bigfloat(-1//8)+polylog2_bigfloat(1//9) ≈ -logBF(9//8)^2/2 atol = ε 
    @test 36*polylog2_bigfloat(1//2)-36*polylog2_bigfloat(1//4)-12*polylog2_bigfloat(1//8)+6*polylog2_bigfloat(1//64) ≈ piBF^2 atol = 32*ε
    @test polylog2_bigfloat(-1) ≈ -piBF^2/12 atol = ε
    @test polylog2_bigfloat(0) == 0.0   
    @test polylog2_bigfloat(1//2) ≈ piBF^2/12 - logBF(2)^2/2 atol = ε
    @test polylog2_bigfloat(1) ≈ piBF^2/6 atol = ε 
    @test polylog2_bigfloat(-1/φBF) ≈ -piBF^2/15 + logBF(φ)^2/2 atol = ε
    @test polylog2_bigfloat(-φBF) ≈ -piBF^2/10 - logBF(φ)^2 atol = ε
    @test polylog2_bigfloat(2-φBF) ≈ piBF^2/15 - logBF(φ)^2 atol = ε
    @test polylog2_bigfloat(1/φBF) ≈ piBF^2/10 - logBF(φ)^2 atol = ε 
    @test polylog2_bigfloat(sqrt(convert(BigFloat,2))-1)-polylog2_bigfloat(1-sqrt(convert(BigFloat,2))) ≈ piBF^2/8 - logBF(1+sqrt(convert(BigFloat,2)))^2/2 atol = 4*ε
    @test polylog2_bigfloat(φBF) ≈ (11*piBF^2)/15 + logBF(Complex(-1/φBF))^2/2 atol = 4*ε
    @test polylog2_bigfloat(φBF^2) ≈ -11*piBF^2/15 - log(Complex(-φBF))^2 atol = 4*ε
end

function spence(x) 
    polylog2_bigfloat(1-x)
end

# Table 27.7 Abramowitz & Stegun. There are four table values with
# a last digit that is off by one digit.

println()
myprintln("Table 27.7 Abramowitz & Stegun")
@testset begin 
   @test spence(0.0)  ≈ 1.644934067 atol=1.0e-9
   @test spence(0.01) ≈ 1.588625448 atol=1.0e-9
   @test spence(0.02) ≈ 1.545799712 atol=1.0e-9
   @test spence(0.03) ≈ 1.507899041 atol=1.0e-9
   @test spence(0.04) ≈ 1.473125860 atol=1.0e-9
   @test spence(0.05) ≈ 1.440633797 atol=1.0e-9
   @test spence(0.06) ≈ 1.409928300 atol=1.0e-9
   @test spence(0.07) ≈ 1.380685041 atol=1.0e-9
   @test spence(0.08) ≈ 1.352675161 atol=1.0e-9
   @test spence(0.09) ≈ 1.325728728 atol=1.0e-9
   @test spence(0.1)  ≈ 1.299714723 atol=1.0e-9
   @test spence(0.11) ≈ 1.274529160 atol=1.0e-9
   @test spence(0.12) ≈ 1.250087584 atol=1.0e-9
   @test spence(0.13) ≈ 1.226320101 atol=1.0e-9
   @test spence(0.14) ≈ 1.203167961 atol=1.0e-9
   @test spence(0.15) ≈ 1.180581124 atol=1.0e-9
   @test spence(0.16) ≈ 1.158516487 atol=1.0e-9 # table wrong last digit
   @test spence(0.17) ≈ 1.136936560 atol=1.0e-9
   @test spence(0.18) ≈ 1.115808451 atol=1.0e-9
   @test spence(0.19) ≈ 1.095103088 atol=1.0e-9
   @test spence(0.20) ≈ 1.074794600 atol=1.0e-9

   @test spence(0.21) ≈ 1.054859830 atol=1.0e-9
   @test spence(0.22) ≈ 1.035277934 atol=1.0e-9
   @test spence(0.23) ≈ 1.016030062 atol=1.0e-9
   @test spence(0.24) ≈ 0.997099088 atol=1.0e-9
   @test spence(0.25) ≈ 0.978469393 atol=1.0e-9
   @test spence(0.26) ≈ 0.960126675 atol=1.0e-9
   @test spence(0.27) ≈ 0.942057798 atol=1.0e-9
   @test spence(0.28) ≈ 0.924250654 atol=1.0e-9
   @test spence(0.29) ≈ 0.906694053 atol=1.0e-9
   @test spence(0.3) ≈ 0.889377624 atol=1.0e-9
   @test spence(0.31) ≈ 0.872291733 atol=1.0e-9
   @test spence(0.32) ≈ 0.855427404 atol=1.0e-9
   @test spence(0.33) ≈ 0.838776261 atol=1.0e-9
   @test spence(0.34) ≈ 0.822330471 atol=1.0e-9
   @test spence(0.35) ≈ 0.806082689 atol=1.0e-9 # table wrong last digit
   @test spence(0.36) ≈ 0.790026024 atol=1.0e-9
   @test spence(0.37) ≈ 0.774153992 atol=1.0e-9
   @test spence(0.38) ≈ 0.758460483 atol=1.0e-9 # table wrong last digit
   @test spence(0.39) ≈ 0.742939737 atol=1.0e-9
   @test spence(0.40) ≈ 0.727586308 atol=1.0e-9
   @test spence(0.41) ≈ 0.712395042 atol=1.0e-9
   @test spence(0.42) ≈ 0.697361058 atol=1.0e-9
   @test spence(0.43) ≈ 0.682479725 atol=1.0e-9
   @test spence(0.44) ≈ 0.667746644 atol=1.0e-9
   @test spence(0.45) ≈ 0.653157631 atol=1.0e-9 # table wrong last digit
   @test spence(0.46) ≈ 0.638708705 atol=1.0e-9
   @test spence(0.47) ≈ 0.624396071 atol=1.0e-9
   @test spence(0.48) ≈ 0.610216108 atol=1.0e-9
   @test spence(0.49) ≈ 0.596165361 atol=1.0e-9
   @test spence(0.50) ≈ 0.582240526 atol=1.0e-9
end

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
 
 # Test dlmf_25_12_3_E3 inside and on unit circle.
 function polylog2_test1(T::DataType,n::Int64)
    OK = true
    for i = 0 : n
        for j = 0 : n
            x = (convert(T, (i/n) * cis(2*pi* j /n)))
            OK = OK && (dlmf_25_12_3_E3(x) < 16)
        end
    end
    OK   
 end

println()
myprintln("Test DLMF identity 25.12.3E3")
@testset begin 
    @test polylog2_test1(Complex{Float16},100) == true
    @test polylog2_test1(Complex{Float32},100) == true
    @test polylog2_test1(Complex{Float64},100) == true
    @test polylog2_test1(Complex{BigFloat},100) == true
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
    OK = true
    for i = 0 : n
        for j = 0 : n
            x = (convert(T, (i/n) * cis(2*pi* j /n)))
            OK = OK && (dlmf_25_12_3_E3(x) < 32)
        end
    end
    OK
 end

println()
myprintln("Test DLMF identity 25.12.E5")
@testset begin 
    @test polylog2_test3(Complex{Float16},100) == true
    @test polylog2_test3(Complex{Float32},100) == true
    @test polylog2_test3(Complex{Float64},100) == true
    @test polylog2_test3(Complex{BigFloat},100) == true
end

#http://dlmf.nist.gov/25.12.E7
function dlmf_25_12_E7(T,x)
    pie = convert(T,pi)
    rd(real(polylog2(cis(x))),  pie*(pie/6 - x/2) + x^2/4)
 end
 
 function polylog2_test4(T::DataType, n::Int64)
    OK = true
    pie = convert(T,pi)
    for i = 0 : n
        x = (2*pie*i)/n
        OK = OK && (dlmf_25_12_E7(T,x) < 32)
    end
    OK
 end

println()
myprintln("Test DLMF identity 25.12.E7")
@testset begin 
    @test polylog2_test4(Float16,100) == true
    @test polylog2_test4(Float32,100) == true
    @test polylog2_test4(Float64,100) == true
    @test polylog2_test4(BigFloat,100) == true
end

#--------------
#See https://en.wikipedia.org/wiki/Dilogarithm
function polylog2_id_1(x)
    rd(polylog2(x)+polylog2(-x),polylog2(x^2)/2) < 16
 end
 
 function polylog2_test5(T::DataType,n::Int64)
       two = convert(T,2)
       OK = true
       for i = 0 : n
        for j = 0 : n
            x = ((two *i)/n) * cis((two*pi*j)/n)
            OK = OK &&  polylog2_id_1(x)
        end
    end
    OK
end

println()
myprintln("Test Dilogarithm reflection identity")
@testset begin 
    @test polylog2_test5(Float16,100) == true
    @test polylog2_test5(Float32,100) == true
    @test polylog2_test5(Float64,100) == true
    @test polylog2_test5(BigFloat,100) == true
end

using PolyLog

function compare_polylog2(T,n)
    OK = true
    tol = 48*eps(T) #seems too big!
    for i = 1 : n
        for j = 0 : n
            x = (convert(T, (i/n) * cis(2*pi* j /n)))
            OK = OK && isapprox(polylog2(x), li2(x), atol=tol, rtol=tol)
        end
    end
  OK      
end


println()
myprintln("Compare to PolyLog")
@testset begin 
    @test compare_polylog2(Complex{Float16},100) == true
    @test compare_polylog2(Complex{Float32},100) == true
    @test compare_polylog2(Complex{Float64},100) == true
    @test compare_polylog2(Complex{BigFloat},100) == true
end

function multiple_prec(x)
    y16 = polylog2(convert(Complex{Float16},x))
    y32 = polylog2(convert(Complex{Float32},x))
    y64 = polylog2(convert(Complex{Float64},x))
    ybf = polylog2(convert(Complex{BigFloat},x))
    isapprox(y16, convert(Complex{Float16},y32), rtol=eps(Float16)) &&
    isapprox(y32, convert(Complex{Float32},y32), rtol=eps(Float32)) &&
    isapprox(y64, convert(Complex{Float64},ybf), rtol=eps(Float64))
end

println()
myprintln("Multiple precision Tests")
@testset begin
    @test multiple_prec(4+8im) == true
    @test multiple_prec(1//4+im//128) == true
    @test multiple_prec(128*im) == true
    @test multiple_prec(1//1 + im//2)
    @test multiple_prec(cis(pi/3)) == true
end

println()   
myprintln("Test Int64 inputs")
@testset begin
    @test polylog2(8) == polylog2(8.0)
    @test polylog2(8 + 5im) == polylog2(8.0 + 5.0im)
end


println()
myprintln("Regression Tests")

@testset begin
  #  polylog2(Float16(-4.53e-6) - Float16(8.3e-7)*im) #1 
  @test polylog2(Float16(-4.53e-6) - Float16(8.3e-7)*im) ==  Float16(-4.53e-6) - Float16(8.3e-7)im
  #  polylog2(0.6 + 0.0im) #3 
  @test  polylog2(0.6 + 0.0im)  == 0.7275863077163334 + 0.0im
  #  polylog2(Float16(0.273) - Float16(0.9004)im) #4 
  @test polylog2(Float16(0.273) - Float16(0.9004)im) == Float16(0.05685) - Float16(0.9385)im
  #  polylog2(2.25 + 0.0im) =/= polylog2(2.25) #5 
  @test polylog2(2.25 + 0.0im) == polylog2(2.25)
end