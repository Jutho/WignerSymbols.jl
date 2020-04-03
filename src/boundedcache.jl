
isnonzero(x::T) where {T<:Integer} = (x != zero(T))
function get_last_nonzero(vec::Vector{T}) where T
    last_nonzero_index = findlast(isnonzero, vec)
    if (isnothing(last_nonzero_index))
        return 1
    end
    return last_nonzero_index
end

# A custom `Integer` subtype to store an integer as its prime factorization
mutable struct StaticPrimeFactorization{U<:Unsigned} <: AbstractPrimeFactorization{U}
    powers::Vector{U}
    sign::Int8
    last_nonzero_index::Int
end
StaticPrimeFactorization(powers::Vector{U}) where {U<:Unsigned} =
    StaticPrimeFactorization{U}(powers, one(Int8), get_last_nonzero(powers))
StaticPrimeFactorization(powers::Vector{U}, sign::Integer) where {U<:Unsigned} =
    StaticPrimeFactorization{U}(powers, sign > 0 ? one(Int8) : -one(Int8), 
        get_last_nonzero(powers))


struct BoundedWignerCache <: AbstractWignerCache
    Wigner3j::Dict{Tuple{UInt,UInt,UInt,Int,Int},Tuple{Rational{BigInt},Rational{BigInt}}}
    Wigner6j::Dict{NTuple{6,UInt},Tuple{Rational{BigInt},Rational{BigInt}}}

    primetable::Array{Int64,1}
    factortable::Array{Array{UInt8,1},1}
    factorialtable::Array{Array{UInt32,1},1}
    bigprimetable::Array{Array{BigInt,1},1}
    bigone::Base.RefValue{BigInt}

    # numerators and denominator buffers in order to GCD the sums
    nums::Array{StaticPrimeFactorization{UInt32},1}
    dens::Array{StaticPrimeFactorization{UInt32},1}

    # prime factorization buffers
    numbuf::StaticPrimeFactorization{UInt32}
    denbuf::StaticPrimeFactorization{UInt32}

    # buffers for r and s calculations. s1n and s1d use numbuf and denbuf.
    s2n::StaticPrimeFactorization{UInt32}  # s2n is reused for sumlist
    snum::StaticPrimeFactorization{UInt32}
    rnum::StaticPrimeFactorization{UInt32}
    sden::StaticPrimeFactorization{UInt32}
    rden::StaticPrimeFactorization{UInt32}
end

factorbuffer(nfact) = StaticPrimeFactorization(zeros(UInt32, nfact), one(Int8),1)


"""
Bounds based on Table 3 of Johansson and Forssén 2016.

Set nj=3 for 3j, nj=6 for 6j, nj=9 for 9j.
"""
function BoundedWignerCache(max_j::T; nj=3) where {T <: Integer}
    multiplier = Int(2 + round(nj / 3))
    maxfactorial = multiplier * max_j + 1
    maxt = max_j + 1

    cache = BoundedWignerCache(
        Dict{Tuple{UInt,UInt,UInt,Int,Int},Tuple{Rational{BigInt},Rational{BigInt}}}(),
        Dict{NTuple{6,UInt},Tuple{Rational{BigInt},Rational{BigInt}}}(),
        [2,3,5],
        [UInt8[], UInt8[1], UInt8[0,1], UInt8[2], UInt8[0,0,1]],
        [UInt32[], UInt32[], UInt32[1], UInt32[1,1], UInt32[3,1], UInt32[3,1,1]],
        [[big(2)], [big(3)], [big(5)]],
        Ref{BigInt}(big(1)),
        [factorbuffer(maxfactorial) for i in 1:maxt],  # nums
        [factorbuffer(maxfactorial) for i in 1:maxt],  # dens
        factorbuffer(maxfactorial),  # numbuf
        factorbuffer(maxfactorial),  # denbuf

        factorbuffer(maxfactorial),  # s2n
        factorbuffer(maxfactorial),  # snum 
        factorbuffer(maxfactorial),  # rnum
        factorbuffer(maxfactorial),  # sden
        factorbuffer(maxfactorial)   # rden
    )
    return cache
end

function precompute_primefactorial!(cache::BoundedWignerCache, n::Integer)
    n < 0 && throw(DomainError(n))
    m = length(cache.factorialtable)-1
    @inbounds while m < n
        prevfactorial = cache.factorialtable[m+1]
        m += 1
        f = primefactor(cache, m).powers
        powers = copy(prevfactorial)
        if length(f) > length(powers) # can at most be 1 larger
            push!(powers, 0)
        end
        for k = 1:length(f)
            powers[k] += f[k]
        end
        push!(cache.factorialtable, powers)
    end
end


"""
Perform an in-place multiplication of a prime factorization of a factorial.
"""
function mul_primefactorial!(cache::BoundedWignerCache, 
        fact_dest::AbstractPrimeFactorization{T}, n::Integer) where T
    precompute_primefactorial!(cache, n)
    row = cache.factorialtable[n+1]
    for i in 1:length(row)
        fact_dest.powers[i] += row[i]
    end
    fact_dest.last_nonzero_index = max(length(row), fact_dest.last_nonzero_index)
    return
end

"""
Convenience function for multiplying together a bunch of primefactorials in a buffer.
"""
function mul_primefactorial!(cache::BoundedWignerCache, 
    fact_dest::AbstractPrimeFactorization{T}, ns::NTuple{N, Integer}) where {T, N}
    for n in ns
        mul_primefactorial!(cache, fact_dest, n)
    end
    return
end

function one!(fact_dest::AbstractPrimeFactorization{T}, sign::Integer) where T
    for i in 1:fact_dest.last_nonzero_index
        fact_dest.powers[i] = zero(T)
    end
    # fill!(fact_dest.powers, zero(T))
    fact_dest.sign = sign < 0 ? -one(Int8) : one(Int8)
    fact_dest.last_nonzero_index = 1
    return
end

# function divgcd!(cache::BoundedWignerCache, 
#         a::AbstractPrimeFactorization, b::AbstractPrimeFactorization)
#     af, bf = a.powers, b.powers
#     last_nonzero_a = findfirst(iszero, af) - 1
#     last_nonzero_b = findfirst(iszero, bf) - 1

#     for k = 1:min(last_nonzero_a, last_nonzero_b)
#         gk = min(af[k], bf[k])
#         af[k] -= gk
#         bf[k] -= gk
#     end
#     return
# end

function _vmin!(a::StaticPrimeFactorization{U}, 
    b::StaticPrimeFactorization{U}) where {U<:Unsigned}
    af = a.powers
    bf = b.powers
    @inbounds for k = 1:a.last_nonzero_index
        af[k] = min(af[k], bf[k])
    end
    return 
end

function _vsub!(a::StaticPrimeFactorization{U}, 
    b::StaticPrimeFactorization{U}) where {U<:Unsigned}
    af = a.powers
    bf = b.powers
    @inbounds for k = 1:b.last_nonzero_index
        af[k] -= bf[k]
    end
    a.last_nonzero_index = max(a.last_nonzero_index, b.last_nonzero_index)

    return a
end

function _vadd!(a::StaticPrimeFactorization{U}, 
    b::StaticPrimeFactorization{U}) where {U<:Unsigned}
    af = a.powers
    bf = b.powers
    @inbounds for k = 1:b.last_nonzero_index
        af[k] = +(af[k], bf[k])
    end
    a.last_nonzero_index = max(a.last_nonzero_index, b.last_nonzero_index)

    return a
end



"""
Given a list of numerators and denominators, compute the common denominator and
the rescaled numerator after putting all fractions at the same common denominator.

This deposits the common denominator in the `cache.denbuf`.
"""
function commondenominator!(cache::BoundedWignerCache, nums::Vector{P}, 
        dens::Vector{P}, nk::Integer) where {P<:AbstractPrimeFactorization}

    isempty(nums) && return one(P)
    iszero(nk) && return one(P)
    # accumulate lcm of denominator

    den = cache.denbuf
    den.powers .= dens[1].powers
    den.last_nonzero_index = dens[1].last_nonzero_index

    for i = 2:nk
        den.powers .= max.(den.powers, dens[i].powers)
        den.last_nonzero_index = max(den.last_nonzero_index, dens[i].last_nonzero_index)
    end

    # rescale numerators
    for i = 1:nk
        _vadd!(nums[i], den)
        _vsub!(nums[i], dens[i])
    end
    return
end

function _convert(cache::BoundedWignerCache, T::Type{BigInt}, a::StaticPrimeFactorization)
    A = one(BigInt)
    # for (n, e) in enumerate(a.powers)
    @inbounds for n in 1:a.last_nonzero_index
        e = a.powers[n]
        if !iszero(e)
            MPZ.mul!(A, bigprime(cache, n, e))
        end
    end
    return a.sign < 0 ? MPZ.neg!(A) : A
end

# auxiliary function to compute sums of a list of PrimeFactorizations as quickly as possible
function sumlist!(cache::BoundedWignerCache, list::Vector{P}, 
        ind = 1:length(list)) where {P <: StaticPrimeFactorization}
    # first compute gcd to take out common factors
    g = cache.snum
    one!(g, 1)
    first_p = list[ind[1]]
    @inbounds for i in 1:first_p.last_nonzero_index
        g.powers[i] = first_p.powers[i]
    end
    g.last_nonzero_index = first_p.last_nonzero_index
    # g = StaticPrimeFactorization(
    #     cache.snum.powers[1:list[ind[1]].last_nonzero_index], one(Int8), 
    #     list[ind[1]].last_nonzero_index)
    # cache.snum.last_nonzero_index = list[ind[1]].last_nonzero_index
    @inbounds for k in ind
        _vmin!(g, list[k])
    end
    @inbounds for k in ind
        _vsub!(list[k], g)
    end
    gint = _convert(cache, BigInt, g)
    L = length(ind)
    if L > 32
        l = L >> 1
        s = sumlist!(cache, list, first(ind).+(0:l-1)) + sumlist!(
                cache, list, first(ind).+(l:L-1))
    else
        # do sum
        s = big(0)
        for k in ind
            MPZ.add!(s, _convert(cache, BigInt, list[k]))
        end
    end
    return MPZ.mul!(s, gint)
end


# squared triangle coefficient
function Δ²!(cache::BoundedWignerCache, s1n::T, s1d::T, 
        j₁, j₂, j₃) where T <: StaticPrimeFactorization
    # also checks the triangle conditions by converting to unsigned integer:
    
    one!(s1n, 1)
    mul_primefactorial!(cache, s1n, (
        convert(UInt, + j₁ + j₂ - j₃),
        convert(UInt, + j₁ - j₂ + j₃),
        convert(UInt, - j₁ + j₂ + j₃))
    )

    one!(s1d, 1)
    mul_primefactorial!(cache, s1d, convert(UInt, j₁ + j₂ + j₃ + 1) )
end

# compute the sum appearing in the 3j symbol
function compute3jseries(cache::BoundedWignerCache, β₁, β₂, β₃, α₁, α₂)
    krange = max(α₁, α₂, zero(α₁)):min(β₁, β₂, β₃)
    numk = length(krange)
    T_int = eltype(eltype(cache.factorialtable))
    T = PrimeFactorization{T_int}

    nums = cache.nums
    dens = cache.dens
    for (i, k) in enumerate(krange)
        num = nums[i]
        den = dens[i]
        one!(num, iseven(k) ? 1 : -1)
        one!(den, 1)
        mul_primefactorial!(cache, den, (k, k-α₁, k-α₂, β₁-k, β₂-k, β₃-k))
        # divgcd!(cache, num, den)
    end

     # write gcd to buffer
    commondenominator!(cache, nums, dens, numk) 
    totalden = _convert(cache, BigInt, cache.denbuf)
    totalnum = sumlist!(cache, nums[1:numk])
    return totalnum//totalden
end

function splitsquare!(s::T, r::T, a::T) where {T <: StaticPrimeFactorization}

    one!(r, a.sign)
    for i in 1:a.last_nonzero_index
        r.powers[i] = convert(UInt8, isodd(a.powers[i]))
    end

    one!(s, 1)
    for i in 1:a.last_nonzero_index
        s.powers[i] = (a.powers[i])>>1
    end

    # maybe this bound could be improved
    r.last_nonzero_index = a.last_nonzero_index
    s.last_nonzero_index = a.last_nonzero_index
end


wigner3j(cache::BoundedWignerCache, j₁, j₂, j₃, m₁, m₂, m₃ = -m₁-m₂) = wigner3j(
    cache, RRBig, j₁, j₂, j₃, m₁, m₂, m₃)
function wigner3j(cache::BoundedWignerCache, T::Type{<:Real}, j₁, j₂, j₃, m₁, m₂, m₃ = -m₁-m₂)
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
        s1n, s1d = cache.numbuf, cache.denbuf
        s2n, snum, rnum, sden, rden = cache.s2n, cache.snum, cache.rnum, cache.sden, cache.rden

        Δ²!(cache, s1n, s1d, j₁, j₂, j₃)
        
        splitsquare!(sden, rden, s1d)
        one!(s2n, 1)
        mul_primefactorial!(cache, s2n, (β₂, β₁ - α₁, β₁ - α₂, β₃, β₃ - α₁, β₂ - α₂))
        s2n_mul_s1n = _vadd!(s1n, s2n)
        splitsquare!(snum, rnum, s2n_mul_s1n)

        
        s = _convert(cache, BigInt, snum) // _convert(cache, BigInt, sden)
        r = _convert(cache, BigInt, rnum) // _convert(cache, BigInt, rden)

        s *= compute3jseries(cache, β₁, β₂, β₃, α₁, α₂)
        cache.Wigner3j[(β₁, β₂, β₃, α₁, α₂)] = (r,s)
    end
    return convert(T, sgn*s)*convert(T, signedroot(r))
end
