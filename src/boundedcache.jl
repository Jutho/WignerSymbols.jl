
struct RSPair{T<:Number} <: Number
    r::T
    s::T
end

const BigNumber = Union{Rational{BigInt},BigInt,BigFloat} 
const Wigner3j = ThreadSafeDict{Tuple{UInt,UInt,UInt,Int,Int},RSPair{Rational{BigInt}}}()
const Wigner6j = ThreadSafeDict{NTuple{6,UInt},RSPair{Rational{BigInt}}}()

function Base.convert(T::Type{<:AbstractFloat}, rs::RSPair)
    convert(T,  # convert to BigFloat first for floating point types
        convert(BigFloat, rs.s)*convert(BigFloat, signedroot(rs.r)))
end
Base.convert(T::Type{<:Real}, rs::RSPair) = convert(T, rs.s)*convert(T, signedroot(rs.r))
function store_result(rs::RSPair{T}) where {T <: BigNumber}
    return deepcopy(rs)
end
store_result(rs::T) where T = rs

function wigner_dicts(Tdict::Type{<:Real})
    return (ThreadSafeDict{Tuple{UInt,UInt,UInt,Int,Int},Tdict}(), 
            ThreadSafeDict{NTuple{6,UInt},Tdict}())
end

function get_thread_caches(Tdict::Type{<:Real}, maxj::Integer)
    d3, d6 = WignerSymbols.wigner_dicts(Float64)
    caches = BoundedWignerCache[]
    resize!(empty!(caches), Threads.nthreads())
    Threads.@threads for i in 1:Threads.nthreads()
        cache = BoundedWignerCache(d3, d6, maxj)
        caches[i] = cache
    end
    return caches
end


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


struct BoundedWignerCache{Tdict} <: AbstractWignerCache

    bounded_multiplier::Int
    max_j::Base.RefValue{Int}

    Wigner3j::ThreadSafeDict{Tuple{UInt,UInt,UInt,Int,Int},Tdict}
    Wigner6j::ThreadSafeDict{NTuple{6,UInt},Tdict}

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

    # buffers for adding
    addbuf::Base.RefValue{BigInt}
    seriesnum::Base.RefValue{BigInt}
    seriesden::Base.RefValue{BigInt}
    
    snumint::Base.RefValue{BigInt}
    sdenint::Base.RefValue{BigInt}
    rnumint::Base.RefValue{BigInt}
    rdenint::Base.RefValue{BigInt}

    # buffers for r and s calculations. s1n and s1d use numbuf and denbuf.
    s2n::StaticPrimeFactorization{UInt32}  # s2n is reused for sumlist
    snum::StaticPrimeFactorization{UInt32}
    rnum::StaticPrimeFactorization{UInt32}
    sden::StaticPrimeFactorization{UInt32}
    rden::StaticPrimeFactorization{UInt32}
end

factorbuffer(nfact) = StaticPrimeFactorization(zeros(UInt32, nfact), one(Int8),1)

"""
Bounds on the prime-counting function.

Rosser, J. Barkley; Schoenfeld, Lowell (1962). "Approximate formulas for some functions of 
prime numbers". Illinois J. Math. 6: 64–94.
"""
πbound(n) = (n < 2) ? 0 : Int(ceil(1.25506n / log(n)))


"""
Stirling's approximation
"""
bound_factorial(n) = Int(ceil(((n+1/2) * log(n) - n ) / log(2)))

"""
Bounds based on Table 3 of Johansson and Forssén 2016.

Set nj=3 for 3j, nj=6 for 6j, nj=9 for 9j.
"""
function BoundedWignerCache(
        dict3j::ThreadSafeDict{Tuple{UInt,UInt,UInt,Int,Int},Tdict}, 
        dict6j::ThreadSafeDict{NTuple{6,UInt},Tdict},
        max_j::Integer, nj=6) where {Tdict}

    multiplier = Int(2 + round(nj / 3))
    maxfactorial = multiplier * max_j + 1
    maxt = max_j + 1
    max_prime_factors = πbound(maxfactorial)

    addbuffer = BigInt(;nbits=bound_factorial(maxfactorial))
    seriesnum = BigInt(;nbits=bound_factorial(maxfactorial))
    seriesden = BigInt(;nbits=bound_factorial(maxfactorial))

    Base.GMP.MPZ.set!(addbuffer, big(0))

    cache = BoundedWignerCache{Tdict}(
        multiplier,
        Ref{Int}(max_j),
        dict3j,
        dict6j,
        [2,3,5],
        [UInt8[], UInt8[1], UInt8[0,1], UInt8[2], UInt8[0,0,1]],
        [UInt32[], UInt32[], UInt32[1], UInt32[1,1], UInt32[3,1], UInt32[3,1,1]],
        [[big(2)], [big(3)], [big(5)]],
        Ref{BigInt}(big(1)),
        [factorbuffer(max_prime_factors) for i in 1:maxt],  # nums
        [factorbuffer(max_prime_factors) for i in 1:maxt],  # dens
        factorbuffer(max_prime_factors),  # numbuf
        factorbuffer(max_prime_factors),  # denbuf

        Ref{BigInt}(addbuffer),
        Ref{BigInt}(seriesnum),
        Ref{BigInt}(seriesden),
        
        Ref{BigInt}(BigInt(;nbits=bound_factorial(maxfactorial))),
        Ref{BigInt}(BigInt(;nbits=bound_factorial(maxfactorial))),
        Ref{BigInt}(BigInt(;nbits=bound_factorial(maxfactorial))),
        Ref{BigInt}(BigInt(;nbits=bound_factorial(maxfactorial))),

        factorbuffer(max_prime_factors),  # s2n
        factorbuffer(max_prime_factors),  # snum 
        factorbuffer(max_prime_factors),  # rnum
        factorbuffer(max_prime_factors),  # sden
        factorbuffer(max_prime_factors)   # rden
    )
    return cache
end

function grow_factorization!(spf::StaticPrimeFactorization{T}, new_size::Integer) where T
    old_size = length(spf.powers)
    resize!(spf.powers, new_size)
    for i in old_size:new_size
        spf.powers[i] = zero(T)
    end
end

function grow!(new_max_j::Integer, cache::BoundedWignerCache)
    cache.max_j[] = new_max_j
    maxfactorial = cache.bounded_multiplier * new_max_j + 1
    maxt = new_max_j + 1
    max_prime_factors = πbound(maxfactorial)

    resize!(cache.nums, maxt)
    resize!(cache.dens, maxt)
    for i in 1:maxt
        cache.nums[i] = factorbuffer(max_prime_factors)
        cache.dens[i] = factorbuffer(max_prime_factors)
    end
    grow_factorization!(cache.numbuf, max_prime_factors)
    grow_factorization!(cache.denbuf, max_prime_factors)
    grow_factorization!(cache.s2n, max_prime_factors)
    grow_factorization!(cache.snum, max_prime_factors)
    grow_factorization!(cache.rnum, max_prime_factors)
    grow_factorization!(cache.sden, max_prime_factors)
    grow_factorization!(cache.rden, max_prime_factors)
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

function _convert!(cache::BoundedWignerCache, A::BigInt, a::StaticPrimeFactorization)
    MPZ.set!(A, cache.bigone[])
    @inbounds for n in 1:a.last_nonzero_index
        e = a.powers[n]
        if !iszero(e)
            MPZ.mul!(A, bigprime(cache, n, e))
        end
    end
    return a.sign < 0 ? MPZ.neg!(A) : A
end

function _convert(cache::BoundedWignerCache, T::Type{BigInt}, a::StaticPrimeFactorization)
    A = one(BigInt)
    return _convert!(cache, A, a)
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
        # summand = big(1)
        for k in ind
            MPZ.add!(s, _convert!(cache, cache.addbuf[], list[k]))
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
    totalden = _convert!(cache, cache.seriesden[], cache.denbuf)
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

# wigner3j(cache::BoundedWignerCache, Wigner3j::Dict, j₁, j₂, j₃, m₁, m₂, m₃ = -m₁-m₂) = wigner3j(
#     cache, RRBig, j₁, j₂, j₃, m₁, m₂, m₃)
function wigner3j(cache::BoundedWignerCache{Tdict}, 
        T::Type{<:Real},
        j₁, j₂, j₃, m₁, m₂, m₃ = -m₁-m₂) where Tdict <: Number
    
    # check angular momenta
    for (jᵢ,mᵢ) in ((j₁, m₁), (j₂, m₂), (j₃, m₃))
        ϵ(jᵢ, mᵢ) || throw(DomainError((jᵢ, mᵢ), "invalid combination (jᵢ, mᵢ)"))
    end
    # check triangle condition and m₁+m₂+m₃ == 0
    if !δ(j₁, j₂, j₃) || !iszero(m₁+m₂+m₃)
        return zero(T)
    end

    new_max_j = Int(ceil(max(abs(j₁), abs(j₂), abs(j₃))))
    if(cache.max_j[] < new_max_j)
        println("GROWING TO $(new_max_j)")
        grow!(2 * new_max_j, cache)
        println("GROWN")
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
        rs = cache.Wigner3j[(β₁, β₂, β₃, α₁, α₂)]
    else
        s1n, s1d = cache.numbuf, cache.denbuf
        s2n, snum, rnum, sden, rden = cache.s2n, cache.snum, cache.rnum, cache.sden, cache.rden

        Δ²!(cache, s1n, s1d, j₁, j₂, j₃)
        splitsquare!(sden, rden, s1d)
        one!(s2n, 1)
        mul_primefactorial!(cache, s2n, (β₂, β₁ - α₁, β₁ - α₂, β₃, β₃ - α₁, β₂ - α₂))
        s2n_mul_s1n = _vadd!(s1n, s2n)
        splitsquare!(snum, rnum, s2n_mul_s1n)

        s = _convert!(cache, cache.snumint[], snum) // _convert!(cache, cache.sdenint[], sden)
        r = _convert!(cache, cache.rnumint[], rnum) // _convert!(cache, cache.rdenint[], rden)

        series = compute3jseries(cache, β₁, β₂, β₃, α₁, α₂)
        # print(typeof(series))
        Base.GMP.MPZ.mul!(s.num, series.num)
        Base.GMP.MPZ.mul!(s.den, series.den)
        rs = convert(Tdict, RSPair(r, s))
        cache.Wigner3j[(β₁, β₂, β₃, α₁, α₂)] = store_result(rs)
    end
    return convert(T, sgn)*convert(T, rs)
end
