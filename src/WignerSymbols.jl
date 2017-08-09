module WignerSymbols
export δ, Δ, clebschgordan, wigner3j, wigner6j

include("primefactorization.jl")

const Wigner3j = Dict{Tuple{UInt,UInt,UInt,Int,Int},Tuple{Rational{BigInt},Rational{BigInt}}}()
const Wigner6j = Dict{NTuple{6,UInt},Tuple{Rational{BigInt},Rational{BigInt}}}()

clebschgordan(j₁, m₁, j₂, m₂, j₃, m₃ = m₁+m₂) = clebschgordan(Float64, j₁, m₁, j₂, m₂, j₃, m₃)
function clebschgordan(T::Type{<:AbstractFloat}, j₁, m₁, j₂, m₂, j₃, m₃ = m₁+m₂)
    s = wigner3j(T, j₁, j₂, j₃, m₁, m₂, -m₃)
    iszero(s) && return s
    if isodd(j₁ - j₂ + m₃)
        s *= -one(s)
    end
    s /= sqrt(convert(T, j₃+j₃+1))
    return s
end

wigner3j(j₁, j₂, j₃, m₁, m₂, m₃ = -m₂-m₃) = wigner3j(Float64, j₁, j₂, j₃, m₁, m₂, m₃)
function wigner3j(T::Type{<:AbstractFloat}, j₁, j₂, j₃, m₁, m₂, m₃ = -m₂-m₃)
    # check angular momenta and triangle condition
    if !(δ(j₁, j₂, j₃) && ϵ(j₁, m₁) && ϵ(j₂, m₂) && ϵ(j₃, m₃))
        throw(DomainError((j₁, j₂, j₃, m₁, m₂, m₃)))
    end
    iszero(m₁+m₂+m₃) || return zero(T)

    # we reorder such that j₁ >= j₂ >= j₃ and m₁ >= 0 or m₁ == 0 && m₂ >= 0
    j₁, j₂, j₃, m₁, m₂, m₃, sgn = reorder3j(j₁, j₂, j₃, m₁, m₂, m₃)
    # do we also want to use Regge symmetries?
    α₁ = convert(Int, j₂ - m₁ - j₃ ) # can be negative
    α₂ = convert(Int, j₁ + m₂ - j₃ ) # can be negative
    β₁ = convert(Int, j₁ + j₂ - j₃ )
    β₂ = convert(Int, j₁ - m₁ )
    β₃ = convert(Int, j₂ + m₂ )

    # extra sign in definition
    sgn = isodd(j₁ - j₃ - m₃) ? -sgn : sgn

    # dictionary lookup or compute
    if haskey(Wigner3j, (β₁, β₂, β₃, α₁, α₂))
        r, s = Wigner3j[(β₁, β₂, β₃, α₁, α₂)]
    else
        s1 = Δ²(j₁, j₂, j₃)
        s2 = prod(map(primefactorial, (β₂, β₁ - α₁, β₁ - α₂, β₃, β₃ - α₁, β₂ - α₂)))

        snum, rnum = splitsquare(s1.num * s2)
        sden, rden = splitsquare(s1.den)
        s = convert(BigInt, snum) // convert(BigInt, sden)
        r = convert(BigInt, rnum) // convert(BigInt, rden)
        s *= compute3jseries(β₁, β₂, β₃, α₁, α₂)
        Wigner3j[(β₁, β₂, β₃, α₁, α₂)] = (r,s)
    end

    return sgn*sqrt(convert(T, r))*convert(T, s)
end

wigner6j(j₁, j₂, j₃, j₄, j₅, j₆) = wigner6j(Float64, j₁, j₂, j₃, j₄, j₅, j₆)
function wigner6j(T::Type{<:AbstractFloat}, j₁, j₂, j₃, j₄, j₅, j₆)
    α̂₁ = (j₁, j₂, j₃)
    α̂₂ = (j₁, j₅, j₆)
    α̂₃ = (j₂, j₄, j₆)
    α̂₄ = (j₃, j₄, j₅)

    # check triangle conditions
    if !(δ(α̂₁...) && δ(α̂₂...) && δ(α̂₃...) && δ(α̂₄...))
        return zero(T)
        # throw(DomainError())
    end
    # reduce
    α₁ = convert(UInt, +(α̂₁...))
    α₂ = convert(UInt, +(α̂₂...))
    α₃ = convert(UInt, +(α̂₃...))
    α₄ = convert(UInt, +(α̂₄...))
    β₁ = convert(UInt, j₁+j₂+j₄+j₅)
    β₂ = convert(UInt, j₁+j₃+j₄+j₆)
    β₃ = convert(UInt, j₂+j₃+j₅+j₆)

    # we should have αᵢ < βⱼ, ∀ i, j and ∑ᵢ αᵢ = ∑ⱼ βⱼ
    # now order them as β₁ >= β₂ >= β₃ >= α₁ >= α₂ >= α₃ >= α₄
    (β₁, β₂, β₃, α₁, α₂, α₃, α₄) = reorder6j(β₁, β₂, β₃, α₁, α₂, α₃, α₄)

    # dictionary lookup or compute
    if haskey(Wigner6j, (β₁, β₂, β₃, α₁, α₂, α₃))
        r, s = Wigner6j[(β₁, β₂, β₃, α₁, α₂, α₃)]
    else
        # order irrelevant: product remains the same under action of reorder6j
        Δ₁ = Δ²(α̂₁...)
        Δ₂ = Δ²(α̂₂...)
        Δ₃ = Δ²(α̂₃...)
        Δ₄ = Δ²(α̂₄...)

        snum, rnum = splitsquare(Δ₁.num * Δ₂.num * Δ₃.num * Δ₄.num)
        sden, rden = splitsquare(Δ₁.den * Δ₂.den * Δ₃.den * Δ₄.den)
        s = convert(BigInt, snum) // convert(BigInt, sden)
        r = convert(BigInt, rnum) // convert(BigInt, rden)
        s *= compute6jseries(β₁, β₂, β₃, α₁, α₂, α₃, α₄)

        Wigner6j[(β₁, β₂, β₃, α₁, α₂, α₃)] = (r, s)
    end

    return sqrt(convert(T, r))*convert(T, s)
end


# check angular momentum
ϵ(j, m) = (abs(m) <= j && isinteger(j-m) && isinteger(j+m))

# check triangle condition
function δ(j₁, j₂, j₃)
    j₃ <= j₁ + j₂ || return false
    j₁ <= j₂ + j₃ || return false
    j₂ <= j₃ + j₁ || return false
    isinteger(j₁ + j₂ - j₃) || return false
    isinteger(j₁ - j₂ + j₃) || return false
    isinteger(j₁ - j₂ - j₃) || return false
    return true
end

# triangle coefficient
function Δ(T::Type{<:AbstractFloat}, j₁, j₂, j₃)
    if !δ(j₁, j₂, j₃)
        throw(DomainError())
    end
    v = Δ²(j₁, j₂, j₃)
    return sqrt(convert(T, convert(BigInt, v.num) // convert(BigInt, v.den)))
end

# squared triangle coefficient
function Δ²(j₁, j₂, j₃)
    # also checks the triangle conditions by converting to unsigned integer:
    n1 = primefactorial( convert(UInt, + j₁ + j₂ - j₃) )
    n2 = primefactorial( convert(UInt, + j₁ - j₂ + j₃) )
    n3 = primefactorial( convert(UInt, - j₁ + j₂ + j₃) )
    d = primefactorial( convert(UInt, j₁ + j₂ + j₃ + 1) )
    # result
    return (n1*n2*n3)//d
end

# reorder parameters determining the 3j symbol to canonical order:
# j₁ >= j₂ >= j₃ and m₁ >= 0 or m₁ == 0 && m₂ >= 0
function reorder3j(j₁, j₂, j₃, m₁, m₂, m₃, sign = one(UInt8))
    if j₁ < j₂
        return reorder3j(j₂, j₁, j₃, m₂, m₁, m₃, -sign)
    elseif j₂ < j₃
        return reorder3j(j₁, j₃, j₂, m₁, m₃, m₂, -sign)
    elseif m₁ < zero(m₁)
        return reorder3j(j₁, j₂, j₃, -m₁, -m₂, -m₃, -sign)
    elseif iszero(m₁) && m₂ < zero(m₂)
        return reorder3j(j₁, j₂, j₃, -m₁, -m₂, -m₃, -sign)
    else
        # sign doesn't matter if total J=j₁ + j₂ + j₃ is even
        if iseven(j₁ + j₂ + j₃)
            sign = one(sign)
        end
        return (j₁, j₂, j₃, m₁, m₂, m₃, sign)
    end
end

# reorder parameters determining the 6j symbol to canonical order:
# β₁ >= β₂ >= β₃ >= α₁ >= α₂ >= α₃ >= α₄
function reorder6j(β₁, β₂, β₃, α₁, α₂, α₃, α₄)
    if β₁ < β₂
        return reorder6j(β₂, β₁, β₃, α₁, α₂, α₃, α₄)
    elseif β₂ < β₃
        return reorder6j(β₁, β₃, β₂, α₁, α₂, α₃, α₄)
    elseif α₁ < α₂
        return reorder6j(β₁, β₂, β₃, α₂, α₁, α₃, α₄)
    elseif α₂ < α₃
        return reorder6j(β₁, β₂, β₃, α₁, α₃, α₂, α₄)
    elseif α₃ < α₄
        return reorder6j(β₁, β₂, β₃, α₁, α₂, α₄, α₃)
    else
        return (β₁, β₂, β₃, α₁, α₂, α₃, α₄)
    end
end

# compute the sum appearing in the 3j symbol
function compute3jseries(β₁, β₂, β₃, α₁, α₂)
    krange = max(α₁, α₂, zero(α₁)):min(β₁, β₂, β₃)
    T = PrimeFactorization{eltype(eltype(factorialtable))}

    nums = Vector{T}(length(krange))
    dens = Vector{T}(length(krange))
    for (i, k) in enumerate(krange)
        num = iseven(k) ? primefactorial(1) : -primefactorial(1)
        den = primefactorial(k)*primefactorial(k-α₁)*primefactorial(k-α₂)*
            primefactorial(β₁-k)*primefactorial(β₂-k)*primefactorial(β₃-k)
        nums[i], dens[i] = divgcd!(num, den)
    end
    den = commondenominator!(nums, dens)
    totalnum = sumlist!(nums)
    totalden = convert(BigInt, PrimeFactorization(den))
    return totalnum//totalden
end

# compute the sum appearing in the 6j symbol
function compute6jseries(β₁, β₂, β₃, α₁, α₂, α₃, α₄)
    krange = max(α₁, α₂, α₃, α₄):min(β₁, β₂, β₃)
    T = PrimeFactorization{eltype(eltype(factorialtable))}

    nums = Vector{T}(length(krange))
    dens = Vector{T}(length(krange))
    for (i, k) in enumerate(krange)
        num = iseven(k) ? primefactorial(k+1) : -primefactorial(k+1)
        den = primefactorial(k-α₁)*primefactorial(k-α₂)*primefactorial(k-α₃)*
            primefactorial(k-α₄)*primefactorial(β₁-k)*primefactorial(β₂-k)*primefactorial(β₃-k)
        nums[i], dens[i] = divgcd!(num, den)
    end
    den = commondenominator!(nums, dens)
    totalnum = sumlist!(nums)
    totalden = convert(BigInt, PrimeFactorization(den))
    return totalnum//totalden
end

end # module
