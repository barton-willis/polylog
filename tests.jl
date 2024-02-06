using Test

φ = Base.MathConstants.golden

function spence(x) 
    polylog2(1-x)
end

# Table 27.7 Abramowitz & Stegun
@testset begin 
    @test spence(0.0) ≈ 1.644934067 atol=1.0e-9
    @test spence(0.01) ≈ 1.588625448 atol=1.0e-9
    @test spence(0.02) ≈ 1.545799712 atol=1.0e-9
    @test spence(0.03) ≈ 1.507899041 atol=1.0e-9
    @test spence(0.04) ≈ 1.473125860 atol=1.0e-9
    @test spence(0.05) ≈ 1.440633797 atol=1.0e-9
    @test spence(0.06) ≈ 1.409928300 atol=1.0e-9
 end
# Special values Int64 See https://en.wikipedia.org/wiki/Dilogarithm
@testset begin  
    @test polylog2(0) == 0
    @test polylog2(1) ≈ pi^2/6
    @test polylog2(2) ≈ pi^2/4 - im*pi*log(2)
end 

function polylog2_f16(x) 
    polylog2(convert(Float16,x))
end
# Special values using binary32 See https://en.wikipedia.org/wiki/Dilogarithm
@testset begin  
    @test polylog2_f16(1/3)-polylog2_f16(1/9)/6 ≈ pi^2/18 - log(3)^2/6
    @test polylog2_f16(-1/3)-polylog2_f16(1/9)/3 ≈ -pi^2/18 + log(3)^2/6
    @test polylog2_f16(-1/2)+polylog2_f16(1/9)/6 ≈ -pi^2/18 + log(2)*log(3)-log(2)^2/2-log(3)^2/3
    @test polylog2_f16(1/4)+polylog2_f16(1/9)/3 ≈ pi^2/18+2*log(2)*log(3)-2*log(2)^2-(2/3)*log(3)^2
    @test polylog2_f16(-1/8)+polylog2_f16(1/9) ≈ -log(9/8)^2/2
    @test 36*polylog2_f16(1/2)-36*polylog2_f16(1/4)-12*polylog2_f16(1/8)+6*polylog2_f16(1/64) ≈ pi^2
    @test polylog2_f16(-1.0) ≈ -pi^2/12
    @test polylog2_f16(0.0) == 0.0   
    @test polylog2_f16(1/2) ≈ pi^2/12 - log(2)^2/2
    @test polylog2_f16(1.0) ≈ pi^2/6
    @test polylog2_f16(-1/φ) ≈ -pi^2/15 + log(φ)^2/2
    @test polylog2_f16(-φ) ≈ -pi^2/10 - log(φ)^2
    @test polylog2_f16(2-φ) ≈ pi^2/15 - log(φ)^2
    @test polylog2_f16(1/φ) ≈ pi^2/10 - log(φ)^2
    @test polylog2_f16(sqrt(2)-1)-polylog2_f16(1-sqrt(2)) ≈ pi^2/8 - log(1+sqrt(2))^2/2 
    @test polylog2_f16(convert(Float64,φ)) ≈ 11*pi^2/15 + log(Complex(-1/φ))^2/2
    @test polylog2_f16(φ^2) ≈ -11*pi^2/15 - log(Complex(-φ))^2
end

function polylog2_f32(x) 
    polylog2(convert(Float32,x))
end

pi32 = convert(Float32,pi)
φ32 = convert(Float32,φ)
# Special values using binary32 See https://en.wikipedia.org/wiki/Dilogarithm
@testset begin  
    @test polylog2_f32(1/3)-polylog2_f32(1/9)/6 ≈ pi32^2/18 - log(3)^2/6 atol=1.0e-7
    @test polylog2_f32(-1/3)-polylog2_f32(1/9)/3 ≈ -pi32^2/18 + log(3)^2/6 atol=1.0e-7
    @test polylog2_f32(-1/2)+polylog2_f32(1/9)/6 ≈ -pi32^2/18 + log(2)*log(3)-log(2)^2/2-log(3)^2/3 atol=1.0e-7
    @test polylog2_f32(1/4)+polylog2_f32(1/9)/3 ≈ pi32^2/18+2*log(2)*log(3)-2*log(2)^2-(2/3)*log(3)^2 atol=1.0e-7
    @test polylog2_f32(-1/8)+polylog2_f32(1/9) ≈ -log(9/8)^2/2 atol=1.0e-8
    @test 36*polylog2_f32(1/2)-36*polylog2_f32(1/4)-12*polylog2_f32(1/8)+6*polylog2_f32(1/64) ≈ pi^2 atol=1.0e-6
    @test polylog2_f32(-1.0) ≈ -pi32^2/12 atol=1.0e-7
    @test polylog2_f32(0.0) == 0.0f0 
    @test polylog2_f32(1/2) ≈ pi32^2/12 - log(convert(Float32,2))^2/2 atol=1.0e-7
    @test polylog2_f32(1.0) ≈ pi32^2/6 atol=1.0e-6
    @test polylog2_f32(-1/φ) ≈ -pi32^2/15 + log(φ32)^2/2 atol=1.0e-7
    @test polylog2_f32(-φ) ≈ -pi32^2/10 - log(φ32)^2 atol=5.0e-7
    @test polylog2_f32(2-φ) ≈ pi32^2/15 - log(φ32)^2 atol=1.0e-7
    @test polylog2_f32(1/φ) ≈ pi32^2/10 - log(φ32)^2 atol=1.0e-7
    @test polylog2_f32(sqrt(2)-1)-polylog2_f32(1-sqrt(2)) ≈ pi32^2/8 - log(1+sqrt(2))^2/2 atol=1.0e-7
    @test polylog2_f32(φ) ≈ 11*pi32^2/15 + clog(-1/φ32)^2/2 atol=5.0e-7
    @test polylog2_f32(φ^2) ≈ -11*pi32^2/15 - clog(-φ32)^2 atol=5.0e-7
end

# Special values using binary64 See https://en.wikipedia.org/wiki/Dilogarithm
@testset begin  
    @test polylog2(1/3)-polylog2(1/9)/6 ≈ pi^2/18 - log(3)^2/6
    @test polylog2(-1/3)-polylog2(1/9)/3 ≈ -pi^2/18 + log(3)^2/6
    @test polylog2(-1/2)+polylog2(1/9)/6 ≈ -pi^2/18 + log(2)*log(3)-log(2)^2/2-log(3)^2/3
    @test polylog2(1/4)+polylog2(1/9)/3 ≈ pi^2/18+2*log(2)*log(3)-2*log(2)^2-(2/3)*log(3)^2
    @test polylog2(-1/8)+polylog2(1/9) ≈ -log(9/8)^2/2
    @test 36*polylog2(1/2)-36*polylog2(1/4)-12*polylog2(1/8)+6*polylog2(1/64) ≈ pi^2
    @test polylog2(-1.0) ≈ -pi^2/12
    @test polylog2(0.0) == 0.0   
    @test polylog2(1/2) ≈ pi^2/12 - log(2)^2/2
    @test polylog2(1.0) ≈ pi^2/6
    @test polylog2(-1/φ) ≈ -pi^2/15 + log(φ)^2/2
    @test polylog2(-φ) ≈ -pi^2/10 - log(φ)^2
    @test polylog2(2-φ) ≈ pi^2/15 - log(φ)^2
    @test polylog2(1/φ) ≈ pi^2/10 - log(φ)^2
    @test polylog2(sqrt(2)-1)-polylog2(1-sqrt(2)) ≈ pi^2/8 - log(1+sqrt(2))^2/2 
    @test polylog2(convert(Float64,φ)) ≈ 11*pi^2/15 + log(Complex(-1/φ))^2/2
    @test polylog2(φ^2) ≈ -11*pi^2/15 - log(Complex(-φ))^2
end

# Table 27.7 Abramowitz & Stegun
@testset begin 
   @test spence(0.0) ≈ 1.644934067 atol=1.0e-9
   @test spence(0.01) ≈ 1.588625448 atol=1.0e-9
   @test spence(0.02) ≈ 1.545799712 atol=1.0e-9
   @test spence(0.03) ≈ 1.507899041 atol=1.0e-9
   @test spence(0.04) ≈ 1.473125860 atol=1.0e-9
   @test spence(0.05) ≈ 1.440633797 atol=1.0e-9
   @test spence(0.06) ≈ 1.409928300 atol=1.0e-9
end