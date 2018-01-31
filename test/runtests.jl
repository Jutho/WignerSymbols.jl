if VERSION < v"0.7.0-DEV.2005"
    const Test = Base.Test
end

using Test
using WignerSymbols
using Compat
using Compat.LinearAlgebra

smalljlist = 0:1//2:10
largejlist = 0:1//2:1000

@testset "triangle coefficient" begin
    for j1 in smalljlist, j2 in smalljlist
        for j3 = abs(j1-j2):(j1+j2)
            @test Δ(j1,j2,j3) ≈ sqrt(factorial(float(j1+j2-j3))*factorial(float(j1-j2+j3))*
                                        factorial(float(j2+j3-j1))/factorial(float(j1+j2+j3+1)))
        end
    end
end

# test 3j:
#--------
@testset "clebschgordan: test orthogonality" begin
    for j1 in smalljlist, j2 in smalljlist
        d1 = 2*j1+1
        d2 = 2*j2+1
        M = zeros(d1*d2, d1*d2)
        ind1 = 1
        for m1 in -j1:j1, m2 in -j2:j2
            ind2 = 1
            for j3 in abs(j1-j2):(j1+j2)
                for m3 in -j3:j3
                    M[ind1,ind2] = clebschgordan(j1,m1,j2,m2,j3,m3)
                    ind2 += 1
                end
            end
            ind1 += 1
        end
        @test M'*M ≈ one(M)
    end
end

# test recurrence relations: Phys Rev E 57, 7274 (1998)
@testset "wigner3j: test recurrence relations" begin
    for k = 1:10
        j2 = convert(BigFloat,rand(0:1//2:1000))
        j3 = convert(BigFloat,rand(0:1//2:1000))
        m2 = convert(BigFloat,rand(-j2:j2))
        m3 = convert(BigFloat,rand(-j3:j3))

        for j in max(abs(j2-j3),abs(m2+m3))+1:(j2+j3)-1
            X = j*sqrt(((j+1)^2-(j2-j3)^2)*((j2+j3+1)^2-(j+1)^2)*((j+1)^2-(m2+m3)^2))
            Y = (2*j+1)*((m2+m3)*(j2*(j2+1)-j3*(j3+1)) - (m2-m3)*j*(j+1))
            Z = (j+1)*sqrt((j^2-(j2-j3)^2)*((j2+j3+1)^2-j^2)*(j^2-(m2+m3)^2))
            tol = 10*max(abs(X),abs(Y),abs(Z))*eps(BigFloat)
            @test (X*wigner3j(BigFloat,j+1,j2,j3,-m2-m3,m2,m3) + Z*wigner3j(BigFloat,j-1,j2,j3,-m2-m3,m2,m3))≈(-Y*wigner3j(BigFloat,j,j2,j3,-m2-m3,m2,m3)) atol=tol
        end
    end
end

# test 6j
#----------
@testset "wigner6j: test orthogonality" begin
    for j1 in smalljlist, j2 in smalljlist, j4 in smalljlist
        for j5 in max(abs(j1-j2-j4),abs(j1-j2+j4),abs(j1+j2-j4)):(j1+j2+j4)
            j6range = max(abs(j2-j4),abs(j1-j5)):min((j2+j4),(j1+j5))
            j3range = max(abs(j1-j2),abs(j4-j5)):min((j1+j2),(j4+j5))
            @test length(j6range) == length(j3range)
            M = zeros(length(j3range),length(j6range))
            for (k2,j6) in enumerate(j6range)
                for (k1,j3) in enumerate(j3range)
                    M[k1,k2] = sqrt(2*j3+1)*sqrt(2*j6+1)*wigner6j(j1,j2,j3,j4,j5,j6)
                end
            end
            @test M'*M ≈ one(M)
        end
    end
end

@testset "wigner6j: test special cases" begin
    for j1 in smalljlist, j2 in smalljlist
        j6 = 0
        j4 = j2
        j5 = j1
        for j3 in abs(j1-j2):(j1+j2)
            @test wigner6j(j1,j2,j3,j4,j5,j6) ≈ (-1)^(j1+j2+j3)/sqrt((2*j1+1)*(2*j2+1))
        end
    end
end

@testset "wigner6j: test recurrence relation" begin
    for k = 1:10
        j2 = convert(BigFloat,rand(largejlist))
        j3 = convert(BigFloat,rand(largejlist))
        l1 = convert(BigFloat,rand(largejlist))
        l2 = convert(BigFloat,rand(abs(l1-j3):(l1+j3)))
        l3 = convert(BigFloat,rand(abs(l1-j2):min(l1+j2)))

        for j in intersect(abs(j2-j3):(j2+j3), abs(l2-l3):(l2+l3))
            X = j*sqrt(((j+1)^2-(j2-j3)^2)*((j2+j3+1)^2-(j+1)^2)*((j+1)^2-(l2-l3)^2)*((l2+l3+1)^2 - (j+1)^2))
            Y = (2*j+1)*( j*(j+1)*( -j*(j+1) + j2*(j2+1) + j3*(j3+1) - 2*l1*(l1+1)) +
                            l2*(l2+1)*( j*(j+1) + j2*(j2+1) - j3*(j3+1) ) +
                            l3*(l3+1)*( j*(j+1) - j2*(j2+1) + j3*(j3+1) ) )
            Z = (j+1)*sqrt((j^2-(j2-j3)^2)*((j2+j3+1)^2-j^2)*(j^2-(l2-l3)^2)*((l2+l3+1)^2 - j^2))
            tol = 10*max(abs(X),abs(Y),abs(Z))*eps(BigFloat)
            @test (X*wigner6j(BigFloat,j+1,j2,j3,l1,l2,l3) + Z*wigner6j(BigFloat,j-1,j2,j3,l1,l2,l3))≈(-Y*wigner6j(BigFloat,j,j2,j3,l1,l2,l3)) atol=tol
        end
    end
end

@testset "test recoupling relation between 3j/clebschgordan and 6j/racahW symbols" begin
    smallerjlist = 0:1//2:5
    for j1 in smallerjlist, j2 in smallerjlist, j3 in smallerjlist
        m1range = -j1:j1
        m2range = -j2:j2
        m3range = -j3:j3
        V1 = Array{Float64}(uninitialized, length(m1range),length(m2range),length(m3range))
        V2 = Array{Float64}(uninitialized, length(m1range),length(m2range),length(m3range))
        for J in max(abs(j1-j2-j3),abs(j1-j2+j3),abs(j1+j2-j3)):(j1+j2+j3)
            J12range = max(abs(j1-j2),abs(J-j3)):min((j1+j2),(J+j3))
            J23range = max(abs(j2-j3),abs(j1-J)):min((j2+j3),(j1+J))
            for J12 in J12range, J23 in J23range
                M = rand(-J:J) # only test for one instance of M in -J:J, should be independent of M anyway
                fill!(V1,0)
                fill!(V2,0)
                for (k1,m1) in enumerate(m1range)
                    for (k2,m2) in enumerate(m2range)
                        abs(m1+m2)<=J12 || continue
                        for (k3,m3) in enumerate(m3range)
                            abs(m2+m3)<=J23 || continue
                            m1+m2+m3==M || continue
                            V1[k1,k2,k3] = clebschgordan(j1,m1,j2,m2,J12)*clebschgordan(J12,m1+m2,j3,m3,J)
                            V2[k1,k2,k3] = clebschgordan(j2,m2,j3,m3,J23)*clebschgordan(j1,m1,J23,m2+m3,J)
                        end
                    end
                end
                @test racahW(j1,j2,J,j3,J12,J23) ≈ vecdot(V2,V1)/sqrt((2*J12+1)*(2*J23+1)) atol=10*eps(Float64)
            end
        end
    end
end
