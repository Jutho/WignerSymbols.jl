using Test
using WignerSymbols
using LinearAlgebra
using Random

Random.seed!(1234)
smalljlist = 0:1//2:10
largejlist = 0:1//2:1000

cache = WignerCache()

# test recurrence relations: Phys Rev E 57, 7274 (1998)
@testset "bounded wigner3j: test recurrence relations" begin
    for k = 1:10
        j2 = convert(BigFloat, rand(0:1//2:1000))
        j3 = convert(BigFloat, rand(0:1//2:1000))
        m2 = convert(BigFloat, rand(-j2:j2))
        m3 = convert(BigFloat, rand(-j3:j3))

        for j in max(abs(j2-j3),abs(m2+m3))+1:(j2+j3)-1
            X = j*sqrt(((j+1)^2-(j2-j3)^2)*((j2+j3+1)^2-(j+1)^2)*((j+1)^2-(m2+m3)^2))
            Y = (2*j+1)*((m2+m3)*(j2*(j2+1)-j3*(j3+1)) - (m2-m3)*j*(j+1))
            Z = (j+1)*sqrt((j^2-(j2-j3)^2)*((j2+j3+1)^2-j^2)*(j^2-(m2+m3)^2))
            tol = 10*max(abs(X),abs(Y),abs(Z))*eps(BigFloat)
            @test (X*wigner3j(cache, BigFloat,j+1,j2,j3,-m2-m3,m2,m3) +
                        Z*wigner3j(cache, BigFloat,j-1,j2,j3,-m2-m3,m2,m3)) ≈
                            (-Y*wigner3j(cache, BigFloat,j,j2,j3,-m2-m3,m2,m3)) atol=tol
        end
    end
end

@testset "bounded wigner3j: test orthogonality relations" begin
    # equivalent to Clebsch-Gordan orthogonality, now test using Float32
    for j1 in smalljlist, j2 in smalljlist
        d1::Int = 2*j1+1
        d2::Int = 2*j2+1
        M = zeros(Float32, (d1*d2, d1*d2))
        ind2 = 1
        for m1 in -j1:j1, m2 in -j2:j2
            ind1 = 1
            @inbounds for j3 in abs(j1-j2):(j1+j2), m3 in -j3:j3
                d3::Int = 2*j3+1
                M[ind1,ind2] += sqrt(d3) * wigner3j(cache, Float32, j1, j2, j3, m1, m2, m3)
                ind1 += 1
            end
            ind2 += 1
        end
        @test M'*M ≈ one(M) # orthogonality relation type 1
        @test M*M' ≈ one(M) # orthogonality relation type 2
    end
end