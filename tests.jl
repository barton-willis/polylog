using Test

# Special values Int64 See https://en.wikipedia.org/wiki/Dilogarithm
@testset begin  
    @test polylog2(0) == 0
    @test polylog2(1) ≈ pi^2/6
    @test polylog2(2) ≈ pi^2/4 - im*pi*log(2)
end 

# Special values binary64 See https://en.wikipedia.org/wiki/Dilogarithm
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
    @test polylog2((1-sqrt(5))/2) ≈ -pi^2/15 + log((1+sqrt(5))/2)^2/2
    @test polylog2(-(1+sqrt(5))/2) ≈ -pi^2/10 - log((1+sqrt(5))/2)^2
    @test polylog2((3-sqrt(5))/2) ≈ pi^2/15 - log((1+sqrt(5))/2)^2
    @test polylog2((sqrt(5)-1)/2) ≈ pi^2/10 - log((1+sqrt(5))/2)^2
end
