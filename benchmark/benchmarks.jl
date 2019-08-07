# some random benchmarks; for j=50000 you need more than 8GB of memory
using WignerSymbols

@time wigner3j(BigFloat, 15, 30, 40, 2, 2, -4)
@time wigner3j(BigFloat, 200, 200, 200, -10, 60, -50)
@time wigner3j(BigFloat, 50000, 50000, 50000, 1000, -6000, 5000)

j = 8
@time wigner6j(BigFloat, j, j, j, j, j, j)
j = 200
@time wigner6j(BigFloat, j, j, j, j, j, j)
j = 600
@time wigner6j(BigFloat, j, j, j, j, j, j)
j = 10000
@time wigner6j(BigFloat, j, j, j, j, j, j)
j = 50000
@time wigner6j(BigFloat, j, j, j, j, j, j)

function compute3jmax(jmax)
    for j₁ = 0:1//2:jmax
        for j₂ = 0:1//2:j₁
            for j₃ = abs(j₁-j₂):j₂
                if isinteger(j₁)
                    m₁ = 0
                    M₂ = min(j₂, j₃-m₁)
                    for m₂ = M₂:(-1):0
                        wigner3j(Float64, j₁,j₂,j₃,0, m₂)
                    end
                    for m₁ = 1:j₁
                        M₂ = min(j₂, j₃-m₁)
                        for m₂ = -M₂:M₂
                            wigner3j(Float64, j₁, j₂, j₃, m₁, m₂)
                        end
                    end
                else
                    for m₁ = 1//2:j₁
                        M₂ = min(j₂, j₃-m₁)
                        for m₂ = -M₂:M₂
                            wigner3j(Float64, j₁, j₂, j₃, m₁, m₂)
                        end
                    end
                end
            end
        end
    end
    return nothing
end

@time compute3jmax(10)

function computeinteger3jmax(jmax)
    for j₁ = 0:jmax
        for j₂ = 0:j₁
            for j₃ = abs(j₁-j₂):j₂
                m₁ = 0
                M₂ = min(j₂, j₃-m₁)
                for m₂ = M₂:(-1):0
                    wigner3j(Float64, j₁,j₂,j₃,0, m₂)
                end
                for m₁ = 1:j₁
                    M₂ = min(j₂, j₃-m₁)
                    for m₂ = -M₂:M₂
                        wigner3j(Float64, j₁, j₂, j₃, m₁, m₂)
                    end
                end
            end
        end
    end
    return nothing
end

function compute6jmax(jmax)
    for j₁ = 0:1//2:jmax
        for j₂ = 0:1//2:j₁
            for j₃ = abs(j₁-j₂):j₂
                for j₄ = 0:1//2:jmax
                    for j₅ = abs(j₃-j₄):min(jmax,j₃+j₄)
                        for j₆ = min(abs(j₂-j₄),abs(j₁-j₅)):min(j₂+j₄,j₁+j₅,jmax)
                            wigner6j(Float64, j₁, j₂, j₃, j₄, j₅, j₆)
                        end
                    end
                end
            end
        end
    end
    return nothing
end
@time compute6jmax(10)


function computeinteger6jmax(jmax)
    for j₁ = 0:jmax
        for j₂ = 0:j₁
            for j₃ = abs(j₁-j₂):j₂
                for j₄ = 0:jmax
                    for j₅ = abs(j₃-j₄):min(jmax,j₃+j₄)
                        for j₆ = min(abs(j₂-j₄),abs(j₁-j₅)):min(j₂+j₄,j₁+j₅,jmax)
                            wigner6j(Float64, j₁, j₂, j₃, j₄, j₅, j₆)
                        end
                    end
                end
            end
        end
    end
    return nothing
end
