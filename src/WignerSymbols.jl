__precompile__(true)
module WignerSymbols
export δ, Δ, clebschgordan, wigner3j, wigner6j, racahV, racahW

if VERSION <= v"0.7.0-DEV.262"
    include("mpz.jl")
    using .MPZ
else
    using Base.GMP.MPZ
end

using Compat

include("halfinteger.jl")
include("primefactorization.jl")

const Wigner3j = Dict{Tuple{UInt,UInt,UInt,Int,Int},Tuple{Rational{BigInt},Rational{BigInt}}}()
const Wigner6j = Dict{NTuple{6,UInt},Tuple{Rational{BigInt},Rational{BigInt}}}()

function __init__()
    global bigone, bigprimetable, Wigner3j, Wigner6j
    bigone[] = big(1)
    bigprimetable[1][1] = big(2)
    bigprimetable[2][1] = big(3)
    bigprimetable[3][1] = big(5)
    Base.rehash!(Wigner3j)
    Base.rehash!(Wigner6j)
end

# check integerness and correctness of (j,m) angular momentum
ϵ(j, m) = (abs(m) <= j && ishalfinteger(j) && isinteger(j-m) && isinteger(j+m))

# check triangle condition
"""
    δ(j₁, j₂, j₃) -> ::Bool

Checks the triangle conditions `j₃ <= j₁ + j₂`, `j₁ <= j₂ + j₃` and `j₂ <= j₃ + j₁`.
"""
function δ(j₁, j₂, j₃)
    j₃ <= j₁ + j₂ || return false
    j₁ <= j₂ + j₃ || return false
    j₂ <= j₃ + j₁ || return false
    isinteger(j₁+j₂+j₃) || return false
    return true
end

# triangle coefficient
"""
    Δ(T::Type{<:AbstractFloat} = Float64, j₁, j₂, j₃) -> ::T

Computes the triangle coefficient
`Δ(j₁, j₂, j₃) = √((j₁+j₂-j₃)!*(j₁-j₂+j₃)!*(j₂+j₃-j₁)! / (j₁+j₂+j₃+1)!)`
as a type `T` floating point number.

Returns `zero(T)` if the triangle condition `δ(j₁, j₂, j₃)` is not satisfied, but
throws a `DomainError` if the `jᵢ`s are not (half)integer
"""
Δ(j₁, j₂, j₃) = Δ(Float64, j₁, j₂, j₃)
function Δ(T::Type{<:AbstractFloat}, j₁, j₂, j₃)
    for jᵢ in (j₁, j₂, j₃)
        (ishalfinteger(jᵢ) && jᵢ >= 0) || throw(DomainError("invalid jᵢ", jᵢ))
    end
    if !δ(j₁, j₂, j₃)
        return zero(T)
    end
    n, d = Δ²(j₁, j₂, j₃)
    return sqrt(convert(T, convert(BigInt, n) // convert(BigInt, d)))
end

"""
    wigner3j(T::Type{<:AbstractFloat} = Float64, j₁, j₂, j₃, m₁, m₂, m₃ = -m₂-m₁) -> ::T

Compute the value of the Wigner-3j symbol
    ⎛ j₁  j₂  j₃ ⎞
    ⎝ m₁  m₂  m₃ ⎠
as a type `T` floating point number. By default, `T = Float64` and `m₃ = -m₁-m₂`.

Returns `zero(T)` if the triangle condition `δ(j₁, j₂, j₃)` is not satisfied, but
throws a `DomainError` if the `jᵢ`s and `mᵢ`s are not (half)integer or `abs(mᵢ) > jᵢ`.
"""
wigner3j(j₁, j₂, j₃, m₁, m₂, m₃ = -m₁-m₂) = wigner3j(Float64, j₁, j₂, j₃, m₁, m₂, m₃)
function wigner3j(T::Type{<:AbstractFloat}, j₁, j₂, j₃, m₁, m₂, m₃ = -m₁-m₂)
    # check angular momenta
    for (jᵢ,mᵢ) in ((j₁, m₁), (j₂, m₂), (j₃, m₃))
        ϵ(jᵢ, mᵢ) || throw(DomainError((jᵢ, mᵢ), "invalid combination (jᵢ, mᵢ)"))
    end
    # check triangle condition and m₁+m₂+m₃ == 0
    if !δ(j₁, j₂, j₃) || !iszero(m₁+m₂+m₃)
        return zero(T)
    end

    # we reorder such that j₁ >= j₂ >= j₃ and m₁ >= 0 or m₁ == 0 && m₂ >= 0
    j₁, j₂, j₃, m₁, m₂, m₃, sgn = reorder3j(j₁, j₂, j₃, m₁, m₂, m₃)
    # TODO: do we also want to use Regge symmetries?
    α₁ = convert(Int, j₂ - m₁ - j₃ ) # can be negative
    α₂ = convert(Int, j₁ + m₂ - j₃ ) # can be negative
    β₁ = convert(Int, j₁ + j₂ - j₃ )
    β₂ = convert(Int, j₁ - m₁ )
    β₃ = convert(Int, j₂ + m₂ )

    # extra sign in definition: α₁ - α₂ = j₁ + m₂ - j₂ + m₁ = j₁ - j₂ + m₃
    sgn = isodd(α₁ - α₂) ? -sgn : sgn

    # dictionary lookup or compute
    if haskey(Wigner3j, (β₁, β₂, β₃, α₁, α₂))
        r, s = Wigner3j[(β₁, β₂, β₃, α₁, α₂)]
    else
        s1n, s1d = Δ²(j₁, j₂, j₃)
        s2n = prod(map(primefactorial, (β₂, β₁ - α₁, β₁ - α₂, β₃, β₃ - α₁, β₂ - α₂)))

        snum, rnum = splitsquare(s1n*s2n)
        sden, rden = splitsquare(s1d)
        s = convert(BigInt, snum) // convert(BigInt, sden)
        r = convert(BigInt, rnum) // convert(BigInt, rden)
        s *= compute3jseries(β₁, β₂, β₃, α₁, α₂)
        Wigner3j[(β₁, β₂, β₃, α₁, α₂)] = (r,s)
    end

    return sgn*sqrt(convert(T, r.num)/convert(T, r.den))*(convert(T, s.num)/convert(T, s.den))
end

"""
    clebschgordan(T::Type{<:AbstractFloat} = Float64, j₁, m₁, j₂, m₂, j₃, m₃ = m₁+m₂) -> ::T

Compute the value of the Clebsch-Gordan coefficient <j₁, m₁; j₂, m₂ | j₃, m₃ >
as a type `T` floating point number. By default, `T = Float64` and `m₃ = m₁+m₂`.

Returns `zero(T)` if the triangle condition `δ(j₁, j₂, j₃)` is not satisfied, but
throws a `DomainError` if the `jᵢ`s and `mᵢ`s are not (half)integer or `abs(mᵢ) > jᵢ`.
"""
clebschgordan(j₁, m₁, j₂, m₂, j₃, m₃ = m₁+m₂) = clebschgordan(Float64, j₁, m₁, j₂, m₂, j₃, m₃)
function clebschgordan(T::Type{<:AbstractFloat}, j₁, m₁, j₂, m₂, j₃, m₃ = m₁+m₂)
    s = wigner3j(T, j₁, j₂, j₃, m₁, m₂, -m₃)
    iszero(s) && return s
    s *= sqrt(convert(T, j₃+j₃+one(j₃)))
    return isodd(convert(Int,j₁ - j₂ + m₃)) ? -s : s
end

"""
    racahV(T::Type{<:AbstractFloat} = Float64, j₁, j₂, j₃, m₁, m₂, m₃ = -m₁-m₂) -> ::T

Compute the value of Racah's V-symbol
V(j₁, j₂, j₃; m₁, m₂, m₃) = (-1)^(-j₁+j₂+j₃) * ⎛ j₁  j₂  j₃ ⎞
                                               ⎝ m₁  m₂  m₃ ⎠
as a type `T` floating point number. By default, `T = Float64` and `m₃ = -m₁-m₂`.

Returns `zero(T)` if the triangle condition `δ(j₁, j₂, j₃)` is not satisfied, but
throws a `DomainError` if the `jᵢ`s and `mᵢ`s are not (half)integer or `abs(mᵢ) > jᵢ`.
"""
racahV(j₁, j₂, j₃, m₁, m₂, m₃ = -m₁-m₂) = racahV(Float64, j₁, j₂, j₃, m₁, m₂, m₃)
function racahV(T::Type{<:AbstractFloat}, j₁, j₂, j₃, m₁, m₂, m₃ = -m₁-m₂)
    s = wigner3j(T, j₁, j₂, j₃, m₁, m₂, m₃)
    return isodd(convert(Int, -j₁ + j₂ + j₃)) ? -s : s
end

"""
    wigner6j(T::Type{<:AbstractFloat} = Float64, j₁, j₂, j₃, j₄, j₅, j₆) -> ::T

Compute the value of the Wigner-6j symbol
    _⎧ j₁  j₂  j₃ ⎫_
     ⎩ j₄  j₅  j₆ ⎭
as a type `T` floating point number. By default, `T = Float64`.

Returns `zero(T)` if any of triangle conditions `δ(j₁, j₂, j₃)`, `δ(j₁, j₆, j₅)`,
`δ(j₂, j₄, j₆)`, `δ(j₃, j₄, j₅)` are not satisfied, but throws a `DomainError` if
the `jᵢ`s are not integer or halfinteger.
"""
wigner6j(j₁, j₂, j₃, j₄, j₅, j₆) = wigner6j(Float64, j₁, j₂, j₃, j₄, j₅, j₆)
function wigner6j(T::Type{<:AbstractFloat}, j₁, j₂, j₃, j₄, j₅, j₆)
    # check validity of `jᵢ`s
    for jᵢ in (j₁, j₂, j₃, j₄, j₅, j₆)
        (ishalfinteger(jᵢ) && jᵢ >= 0) || throw(DomainError("invalid jᵢ", jᵢ))
    end

    α̂₁ = map(converthalfinteger, (j₁, j₂, j₃))
    α̂₂ = map(converthalfinteger, (j₁, j₆, j₅))
    α̂₃ = map(converthalfinteger, (j₂, j₄, j₆))
    α̂₄ = map(converthalfinteger, (j₃, j₄, j₅))

    # check triangle conditions
    if !(δ(α̂₁...) && δ(α̂₂...) && δ(α̂₃...) && δ(α̂₄...))
        return zero(T)
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
        n₁, d₁ = Δ²(α̂₁...)
        n₂, d₂ = Δ²(α̂₂...)
        n₃, d₃ = Δ²(α̂₃...)
        n₄, d₄ = Δ²(α̂₄...)

        snum, rnum = splitsquare(n₁ * n₂ * n₃ * n₄)
        sden, rden = splitsquare(d₁ * d₂ * d₃ * d₄)
        s = convert(BigInt, snum) // convert(BigInt, sden)
        r = convert(BigInt, rnum) // convert(BigInt, rden)
        s *= compute6jseries(β₁, β₂, β₃, α₁, α₂, α₃, α₄)

        Wigner6j[(β₁, β₂, β₃, α₁, α₂, α₃)] = (r, s)
    end

    return sqrt(convert(T, r.num)/convert(T, r.den))*(convert(T, s.num)/convert(T, s.den))
end

"""
    racahW(T::Type{<:AbstractFloat} = Float64, j₁, j₂, J, j₃, J₁₂, J₂₃) -> ::T

Compute the value of Racah's W coefficient
`W(j₁, j₂, J, j₃; J₁₂, J₂₃) = <(j₁,(j₂j₃)J₂₃)J | ((j₁j₂)J₁₂,j₃)J> / sqrt((2J₁₂+1)*(2J₁₃+1))`
as a type `T` floating point number. By default, `T = Float64`.

Returns `zero(T)` if any of triangle conditions `δ(j₁, j₂, J₁₂)`, `δ(j₂, j₃, J₂₃)`,
`δ(j₁, J₂₃, J)`, `δ(J₁₂, j₃, J)` are not satisfied, but throws a `DomainError` if
the `jᵢ`s and `J`s are not integer or halfinteger.
"""
racahW(j₁, j₂, J, j₃, J₁₂, J₂₃) = racahW(Float64, j₁, j₂, J, j₃, J₁₂, J₂₃)
function racahW(T::Type{<:AbstractFloat}, j₁, j₂, J, j₃, J₁₂, J₂₃)
    s = wigner6j(T, j₁, j₂, J₁₂, j₃, J, J₂₃)
    if !iszero(s) && isodd(convert(Int, j₁ + j₂ + j₃ + J))
        return -s
    else
        return s
    end
end

# COMPUTATIONAL ROUTINES
#------------------------
# squared triangle coefficient
function Δ²(j₁, j₂, j₃)
    # also checks the triangle conditions by converting to unsigned integer:
    n1 = primefactorial( convert(UInt, + j₁ + j₂ - j₃) )
    n2 = primefactorial( convert(UInt, + j₁ - j₂ + j₃) )
    n3 = primefactorial( convert(UInt, - j₁ + j₂ + j₃) )
    d = primefactorial( convert(UInt, j₁ + j₂ + j₃ + 1) )
    # result
    return (n1*n2*n3), d
end

# reorder parameters determining the 3j symbol to canonical order:
# j₁ >= j₂ >= j₃ and m₁ >= 0 or m₁ == 0 && m₂ >= 0
function reorder3j(j₁, j₂, j₃, m₁, m₂, m₃, sign = one(Int8))
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
        if iseven(convert(UInt,j₁ + j₂ + j₃))
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

    nums = Vector{T}(uninitialized, length(krange))
    dens = Vector{T}(uninitialized, length(krange))
    for (i, k) in enumerate(krange)
        num = iseven(k) ? one(T) : -one(T)
        den = primefactorial(k)*primefactorial(k-α₁)*primefactorial(k-α₂)*
            primefactorial(β₁-k)*primefactorial(β₂-k)*primefactorial(β₃-k)
        nums[i], dens[i] = divgcd!(num, den)
    end
    den = commondenominator!(nums, dens)
    totalnum = sumlist!(nums)
    totalden = convert(BigInt, den)
    return totalnum//totalden
end

# compute the sum appearing in the 6j symbol
function compute6jseries(β₁, β₂, β₃, α₁, α₂, α₃, α₄)
    krange = max(α₁, α₂, α₃, α₄):min(β₁, β₂, β₃)
    T = PrimeFactorization{eltype(eltype(factorialtable))}

    nums = Vector{T}(uninitialized, length(krange))
    dens = Vector{T}(uninitialized, length(krange))
    for (i, k) in enumerate(krange)
        num = iseven(k) ? primefactorial(k+1) : -primefactorial(k+1)
        den = primefactorial(k-α₁)*primefactorial(k-α₂)*primefactorial(k-α₃)*
            primefactorial(k-α₄)*primefactorial(β₁-k)*primefactorial(β₂-k)*primefactorial(β₃-k)
        nums[i], dens[i] = divgcd!(num, den)
    end
    den = commondenominator!(nums, dens)
    totalnum = sumlist!(nums)
    totalden = convert(BigInt, den)
    return totalnum//totalden
end

end # module
