__precompile__(true)
module WignerSymbols
export δ, Δ, clebschgordan, wigner3j, wigner6j, racahV, racahW, HalfInteger
export WignerCache, BoundedWignerCache
using Base.GMP.MPZ
using HalfIntegers
using RationalRoots
const RRBig = RationalRoot{BigInt}
import RationalRoots: _convert

# debug stuff
using TimerOutputs
const to = TimerOutput()

include("wignercache.jl")
include("primefactorization.jl")
include("boundedcache.jl")

# imitating Base.Random.THREAD_RNGs, initialize new caches on each thread
const wigner_caches = WignerCache[]
@inline get_local_cache() = get_local_cache(Threads.threadid())
@noinline function get_local_cache(tid::Int)
    @assert 0 < tid <= length(wigner_caches)
    if @inbounds isassigned(wigner_caches, tid)
        @inbounds cache = wigner_caches[tid]
    else
        cache = WignerCache()
        @inbounds wigner_caches[tid] = cache
    end
    return cache
end

function __init__()
    global wigner_caches
    resize!(empty!(wigner_caches), Threads.nthreads())
end

# check integerness and correctness of (j,m) angular momentum
ϵ(j, m) = (abs(m) <= j && ishalfinteger(j) && isinteger(j-m) && isinteger(j+m))

# check triangle condition
"""
    δ(j₁, j₂, j₃) -> ::Bool

Checks the triangle conditions `j₃ <= j₁ + j₂`, `j₁ <= j₂ + j₃` and `j₂ <= j₃ + j₁`.
"""
δ(j₁, j₂, j₃) = (j₃ <= j₁ + j₂) && (j₁ <= j₂ + j₃) && (j₂ <= j₃ + j₁) && isinteger(j₁+j₂+j₃)

# triangle coefficient
"""
    Δ(T::Type{<:Real} = RationalRoot{BigInt}, j₁, j₂, j₃) -> ::T

Computes the triangle coefficient
`Δ(j₁, j₂, j₃) = √((j₁+j₂-j₃)!*(j₁-j₂+j₃)!*(j₂+j₃-j₁)! / (j₁+j₂+j₃+1)!)`
as a type `T` real number.

Returns `zero(T)` if the triangle condition `δ(j₁, j₂, j₃)` is not satisfied, but
throws a `DomainError` if the `jᵢ`s are not (half)integer
"""
Δ(j₁, j₂, j₃) = Δ(get_local_cache(), RRBig, j₁, j₂, j₃)
Δ(cache::WignerCache, j₁, j₂, j₃) = Δ(cache, RRBig, j₁, j₂, j₃)
function Δ(cache::WignerCache, T::Type{<:Real}, j₁, j₂, j₃)
    for jᵢ in (j₁, j₂, j₃)
        (ishalfinteger(jᵢ) && jᵢ >= zero(jᵢ)) || throw(DomainError("invalid jᵢ", jᵢ))
    end
    if !δ(j₁, j₂, j₃)
        return zero(T)
    end
    n, d = Δ²(cache, j₁, j₂, j₃)
    return convert(T,
        signedroot(RationalRoot{BigInt},
            _convert(cache, BigInt, n)//_convert(cache, BigInt, d)))
end

"""
    wigner3j(T::Type{<:Real} = RationalRoot{BigInt}, j₁, j₂, j₃, m₁, m₂, m₃ = -m₂-m₁) -> ::T

Compute the value of the Wigner-3j symbol
    ⎛ j₁  j₂  j₃ ⎞
    ⎝ m₁  m₂  m₃ ⎠
as a type `T` real number. By default, `T = RationalRoot{BigInt}` and `m₃ = -m₁-m₂`.

Returns `zero(T)` if the triangle condition `δ(j₁, j₂, j₃)` is not satisfied, but
throws a `DomainError` if the `jᵢ`s and `mᵢ`s are not (half)integer or `abs(mᵢ) > jᵢ`.
"""
wigner3j(j₁::Integer, j₂, j₃, m₁, m₂, m₃ = -m₁-m₂) = wigner3j(
    get_local_cache(), RRBig, j₁, j₂, j₃, m₁, m₂, m₃)
wigner3j(T::Type{<:Real}, j₁, j₂, j₃, m₁, m₂, m₃ = -m₁-m₂) = wigner3j(
    get_local_cache(), T, j₁, j₂, j₃, m₁, m₂, m₃)
wigner3j(cache::WignerCache, j₁, j₂, j₃, m₁, m₂, m₃ = -m₁-m₂) = wigner3j(
    cache, RRBig, j₁, j₂, j₃, m₁, m₂, m₃)
function wigner3j(cache::WignerCache, T::Type{<:Real}, j₁, j₂, j₃, m₁, m₂, m₃ = -m₁-m₂)
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
    if haskey(cache.Wigner3j, (β₁, β₂, β₃, α₁, α₂))
        r, s = cache.Wigner3j[(β₁, β₂, β₃, α₁, α₂)]
    else
        s1n, s1d = Δ²(cache, j₁, j₂, j₃)
        s2n = (
            primefactorial(cache, β₂) * primefactorial(cache, β₁ - α₁) * primefactorial(cache, β₁ - α₂) *
            primefactorial(cache, β₃) * primefactorial(cache, β₃ - α₁) * primefactorial(cache, β₂ - α₂))

        snum, rnum = splitsquare(s1n*s2n)
        sden, rden = splitsquare(s1d)

        s = _convert(cache, BigInt, snum) // _convert(cache, BigInt, sden)
        r = _convert(cache, BigInt, rnum) // _convert(cache, BigInt, rden)
        s *= compute3jseries(cache, β₁, β₂, β₃, α₁, α₂)
        cache.Wigner3j[(β₁, β₂, β₃, α₁, α₂)] = (r,s)
    end
    return convert(T, sgn*s)*convert(T, signedroot(r))
end

"""
    clebschgordan(T::Type{<:Real} = RationalRoot{BigInt}, j₁, m₁, j₂, m₂, j₃, m₃ = m₁+m₂) -> ::T

Compute the value of the Clebsch-Gordan coefficient <j₁, m₁; j₂, m₂ | j₃, m₃ >
as a type `T` real number. By default, `T = RationalRoot{BigInt}` and `m₃ = m₁+m₂`.

Returns `zero(T)` if the triangle condition `δ(j₁, j₂, j₃)` is not satisfied, but
throws a `DomainError` if the `jᵢ`s and `mᵢ`s are not (half)integer or `abs(mᵢ) > jᵢ`.
"""
clebschgordan(j₁, m₁, j₂, m₂, j₃, m₃ = m₁+m₂) =
    clebschgordan(get_local_cache(), RRBig, j₁, m₁, j₂, m₂, j₃, m₃)
clebschgordan(cache::WignerCache, j₁, m₁, j₂, m₂, j₃, m₃ = m₁+m₂) =
    clebschgordan(cache, RRBig, j₁, m₁, j₂, m₂, j₃, m₃)
function clebschgordan(cache::WignerCache, T::Type{<:Real}, j₁, m₁, j₂, m₂, j₃, m₃ = m₁+m₂)
    s = wigner3j(cache, T, j₁, j₂, j₃, m₁, m₂, -m₃)
    iszero(s) && return s
    s *= convert(T, signedroot(RRBig, j₃+j₃+one(j₃)))
    return isodd(convert(Int,j₁ - j₂ + m₃)) ? -s : s
end

"""
    racahV(T::Type{<:Real} = RationalRoot{BigInt}, j₁, j₂, j₃, m₁, m₂, m₃ = -m₁-m₂) -> ::T

Compute the value of Racah's V-symbol
V(j₁, j₂, j₃; m₁, m₂, m₃) = (-1)^(-j₁+j₂+j₃) * ⎛ j₁  j₂  j₃ ⎞
                                               ⎝ m₁  m₂  m₃ ⎠
as a type `T` real number. By default, `T = RationalRoot{BigInt}` and `m₃ = -m₁-m₂`.

Returns `zero(T)` if the triangle condition `δ(j₁, j₂, j₃)` is not satisfied, but
throws a `DomainError` if the `jᵢ`s and `mᵢ`s are not (half)integer or `abs(mᵢ) > jᵢ`.
"""
racahV(j₁, j₂, j₃, m₁, m₂, m₃ = -m₁-m₂) = racahV(
    get_local_cache(), RRBig, j₁, j₂, j₃, m₁, m₂, m₃)
racahV(cache::WignerCache, j₁, j₂, j₃, m₁, m₂, m₃ = -m₁-m₂) = racahV(
    cache, RRBig, j₁, j₂, j₃, m₁, m₂, m₃)
function racahV(cache::WignerCache, T::Type{<:Real}, j₁, j₂, j₃, m₁, m₂, m₃ = -m₁-m₂)
    s = wigner3j(cache, T, j₁, j₂, j₃, m₁, m₂, m₃)
    return isodd(convert(Int, -j₁ + j₂ + j₃)) ? -s : s
end

"""
    wigner6j(T::Type{<:Real} = RationalRoot{BigInt}, j₁, j₂, j₃, j₄, j₅, j₆) -> ::T

Compute the value of the Wigner-6j symbol
    _⎧ j₁  j₂  j₃ ⎫_
     ⎩ j₄  j₅  j₆ ⎭
as a type `T` real number. By default, `T = RationalRoot{BigInt}`.

Returns `zero(T)` if any of triangle conditions `δ(j₁, j₂, j₃)`, `δ(j₁, j₆, j₅)`,
`δ(j₂, j₄, j₆)`, `δ(j₃, j₄, j₅)` are not satisfied, but throws a `DomainError` if
the `jᵢ`s are not integer or halfinteger.
"""
wigner6j(j₁, j₂, j₃, j₄, j₅, j₆) = wigner6j(get_local_cache(), RRBig, j₁, j₂, j₃, j₄, j₅, j₆)
wigner6j(T::Type{<:Real}, j₁, j₂, j₃, j₄, j₅, j₆) = wigner6j(get_local_cache(), T, j₁, j₂, j₃, j₄, j₅, j₆)
wigner6j(cache::WignerCache, j₁, j₂, j₃, j₄, j₅, j₆) = wigner6j(cache, RRBig, j₁, j₂, j₃, j₄, j₅, j₆)
function wigner6j(cache::WignerCache, T::Type{<:Real}, j₁, j₂, j₃, j₄, j₅, j₆)
    # check validity of `jᵢ`s
    for jᵢ in (j₁, j₂, j₃, j₄, j₅, j₆)
        (ishalfinteger(jᵢ) && jᵢ >= zero(jᵢ)) || throw(DomainError("invalid jᵢ", jᵢ))
    end

    α̂₁ = (j₁, j₂, j₃)
    α̂₂ = (j₁, j₆, j₅)
    α̂₃ = (j₂, j₄, j₆)
    α̂₄ = (j₃, j₄, j₅)

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
    if haskey(cache.Wigner6j, (β₁, β₂, β₃, α₁, α₂, α₃))
        r, s = cache.Wigner6j[(β₁, β₂, β₃, α₁, α₂, α₃)]
    else
        # order irrelevant: product remains the same under action of reorder6j
        n₁, d₁ = Δ²(cache, α̂₁...)
        n₂, d₂ = Δ²(cache, α̂₂...)
        n₃, d₃ = Δ²(cache, α̂₃...)
        n₄, d₄ = Δ²(cache, α̂₄...)

        snum, rnum = splitsquare(n₁ * n₂ * n₃ * n₄)
        sden, rden = splitsquare(d₁ * d₂ * d₃ * d₄)
        snu, sden = divgcd!(snum, sden)
        rnu, rden = divgcd!(rnum, rden)
        s = _convert(cache, BigInt, snum) // _convert(cache, BigInt, sden)
        r = _convert(cache, BigInt, rnum) // _convert(cache, BigInt, rden)
        s *= compute6jseries(cache, β₁, β₂, β₃, α₁, α₂, α₃, α₄)

        cache.Wigner6j[(β₁, β₂, β₃, α₁, α₂, α₃)] = (r, s)
    end
    return _convert(T, s) * convert(T, signedroot(r))
end

"""
    racahW(T::Type{<:Real} = RationalRoot{BigInt}, j₁, j₂, J, j₃, J₁₂, J₂₃) -> ::T

Compute the value of Racah's W coefficient
`W(j₁, j₂, J, j₃; J₁₂, J₂₃) = <(j₁,(j₂j₃)J₂₃)J | ((j₁j₂)J₁₂,j₃)J> / sqrt((2J₁₂+1)*(2J₁₃+1))`
as a type `T` real number. By default, `T = RationalRoot{BigInt}`.

Returns `zero(T)` if any of triangle conditions `δ(j₁, j₂, J₁₂)`, `δ(j₂, j₃, J₂₃)`,
`δ(j₁, J₂₃, J)`, `δ(J₁₂, j₃, J)` are not satisfied, but throws a `DomainError` if
the `jᵢ`s and `J`s are not integer or halfinteger.
"""
racahW(j₁, j₂, J, j₃, J₁₂, J₂₃) = racahW(
    get_local_cache(), RRBig, j₁, j₂, J, j₃, J₁₂, J₂₃)
racahW(cache::WignerCache, j₁, j₂, J, j₃, J₁₂, J₂₃) = racahW(
    cache, RRBig, j₁, j₂, J, j₃, J₁₂, J₂₃)
function racahW(cache::WignerCache, T::Type{<:Real}, j₁, j₂, J, j₃, J₁₂, J₂₃)
    s = wigner6j(cache, T, j₁, j₂, J₁₂, j₃, J, J₂₃)
    if !iszero(s) && isodd(convert(Int, j₁ + j₂ + j₃ + J))
        return -s
    else
        return s
    end
end

# COMPUTATIONAL ROUTINES
#------------------------
# squared triangle coefficient
function Δ²(cache::AbstractWignerCache, j₁, j₂, j₃)
    # also checks the triangle conditions by converting to unsigned integer:
    n1 = primefactorial(cache, convert(UInt, + j₁ + j₂ - j₃) )
    n2 = primefactorial(cache, convert(UInt, + j₁ - j₂ + j₃) )
    n3 = primefactorial(cache, convert(UInt, - j₁ + j₂ + j₃) )
    d = primefactorial(cache, convert(UInt, j₁ + j₂ + j₃ + 1) )
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
function compute3jseries(cache::WignerCache, β₁, β₂, β₃, α₁, α₂)
    krange = max(α₁, α₂, zero(α₁)):min(β₁, β₂, β₃)
    T = PrimeFactorization{eltype(eltype(cache.factorialtable))}

    nums = Vector{T}(undef, length(krange))
    dens = Vector{T}(undef, length(krange))
    @timeit to "3jseries k loop" for (i, k) in enumerate(krange)
        num = iseven(k) ? one(T) : -one(T)
        den = primefactorial(cache, k)*primefactorial(cache, k-α₁)*
            primefactorial(cache, k-α₂)*primefactorial(cache, β₁-k)*
            primefactorial(cache, β₂-k)*primefactorial(cache, β₃-k)
        nums[i], dens[i] = divgcd!(num, den)
    end
    @timeit to "commondenominator!" den = commondenominator!(nums, dens)
    @timeit to "totalnum" totalnum = sumlist!(cache, nums)
    @timeit to "totalden" totalden = _convert(cache, BigInt, den)
    return totalnum//totalden
end

# compute the sum appearing in the 6j symbol
function compute6jseries(cache::WignerCache, β₁, β₂, β₃, α₁, α₂, α₃, α₄)
    krange = max(α₁, α₂, α₃, α₄):min(β₁, β₂, β₃)
    T = PrimeFactorization{eltype(eltype(cache.factorialtable))}

    nums = Vector{T}(undef, length(krange))
    dens = Vector{T}(undef, length(krange))
    for (i, k) in enumerate(krange)
        num = iseven(k) ? primefactorial(cache, k+1) : -primefactorial(cache, k+1)
        den = primefactorial(cache, k-α₁)*primefactorial(cache, k-α₂)*
            primefactorial(cache, k-α₃)*primefactorial(cache, k-α₄)*
            primefactorial(cache, β₁-k)*primefactorial(cache, β₂-k)*primefactorial(cache, β₃-k)
        nums[i], dens[i] = divgcd!(num, den)
    end
    den = commondenominator!(nums, dens)
    totalnum = sumlist!(cache, nums)
    totalden = _convert(cache, BigInt, den)
    return totalnum//totalden
end

end # module
