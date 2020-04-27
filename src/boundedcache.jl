"""
This part of the code primarily dispatches on the `BoundedWignerCache` type on the first
argument of the internal functions. It reduces the number of allocations in heavy workloads
by more than an order of magnitude with a number of different optimizations, compared to the
code in main code file. This makes it accelerate computations with parallelism like
threading without being limited by the memory bandwidth.

* The number n! has a prime factorization containing all prime numbers less than n. Then
  the total number of prime factors one needs to store is bounded by the prime counting 
  function π(n). We thus preallocate prime factorization buffers of size ~ π(n).
* Instead of modifying array sizes with push/pop, one keeps track of the largest nonzero 
  prime factor during additions. Prime factorizations and sums are preallocated with bounds
  given in Johansson and Forssén 2016.
* BigInt is basically an array allocated on the heap, and thus generates a lot of
  painful allocations. Each thread reuses the same BigInts for different parts of the 
  calculation, when possible.
"""

"""
The nj symbols can be written in the form r ⋅ √s, and we store such numbers in the result 
caches `Wigner3j` and `Wigner6j`. Previously this was a `Tuple{T,T}`, but this was 
challenging to use as a Number type. The default result key is `RSPair{Rational{BigInt}}``.
"""
struct RSPair{T<:Number} <: Number
    r::T
    s::T
end

"""
Convert a number of the form r ⋅ √s to some kind of floating point number. This needs to
be done by way of BigFloat to avoid overflows.
"""
function Base.convert(T::Type{<:AbstractFloat}, rs::RSPair)
    convert(T,  # convert to BigFloat first for floating point types
        convert(BigFloat, rs.s)*convert(BigFloat, signedroot(rs.r)))
end

"""
Convert a number of the form r ⋅ √s to an arbitrary number type involves just performing
the conversion and then the actual computation. 
"""
Base.convert(T::Type{<:Real}, rs::RSPair) = convert(T, rs.s)*convert(T, signedroot(rs.r))

"""
The arbitrary-precision types in Julia are still unpolished, and they are unlike normal
numerical types as they mutable. This BigNumber type allows us to dispatch based on 
mutability and make a deepcopy if it's a mutable type.
"""
const BigNumber = Union{Rational{BigInt},BigInt,BigFloat} 

"""
We're trying to reuse BigInts when possible, but this means that if one wants to keep the 
exact solutions in the dictionary, one must deepcopy the number.
"""
store_result(rs::RSPair{T}) where {T <: BigNumber} = deepcopy(rs)
store_result(rs) = rs

"""
The default global thread-safe dictionaries used to cache results of calculations.
They map the symmetrized form to r,s pairs.
"""
const Wigner3j = ThreadSafeDict{Tuple{UInt,UInt,UInt,Int,Int},RSPair{Rational{BigInt}}}()
const Wigner6j = ThreadSafeDict{NTuple{6,UInt},RSPair{Rational{BigInt}}}()

"""
Convenience function for generating thread-safe dictionaries which store the results.
Returns the dictionaries for Wigner3j and Wigner6j in a tuple.
"""
function wigner_dicts(Tdict::Type{<:Real})
    return (ThreadSafeDict{Tuple{UInt,UInt,UInt,Int,Int},Tdict}(),  # 3j
            ThreadSafeDict{NTuple{6,UInt},Tdict}())                 # 6j
end

"""
Convenience function for generating an array of caches for each thread. The default 
behavior essentially calls `get_thread_caches(Rational{BigInt}, MAX_J[])` for the default
caches. 
"""
function get_thread_caches(Tdict::Type{<:Real}, maxj::Integer)
    d3, d6 = WignerSymbols.wigner_dicts(Tdict)  # shared by all threads
    caches = BoundedWignerCache[]
    resize!(empty!(caches), Threads.nthreads())
    Threads.@threads for i in 1:Threads.nthreads()
        cache = BoundedWignerCache(d3, d6, maxj)
        caches[i] = cache
    end
    return caches
end

# utilities for StaticPrimeFactorization
isnonzero(x::T) where {T<:Integer} = (x != zero(T))
function get_last_nonzero(vec::Vector{T}) where T
    last_nonzero_index = findlast(isnonzero, vec)
    if (isnothing(last_nonzero_index))
        return 1
    end
    return last_nonzero_index
end

"""
A custom `Integer` subtype to store an integer as its prime factorization.

Unlike PrimeFactorization, this type keeps track of the index of the last nonzero
index, and it is not expected for the `powers` vector to change size. If more powers are
required, we reallocate every buffer (see `grow!`). When one computes with this number,
all loops over the factors stop at the last nonzero index, which buys ~20% in performance.
"""
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


"""
This structure keeps the buffers into which prime factorizations, big integers, and Wigner
symbols are buffered. By allocating these buffers, one avoids having to reallocate for
each combination of jᵢ.

We parametrize on the type of the stored result `Tdict`.
"""
struct BoundedWignerCache{Tdict} <: AbstractWignerCache

    bounded_multiplier::Int  # 3 for 3j, 4 for 6j. determines the largets factorial
    max_j::Base.RefValue{Int}  # largest jᵢ allowed by the cache

    # the caches for each thread should all point to the same ThreadSafeDict
    Wigner3j::ThreadSafeDict{Tuple{UInt,UInt,UInt,Int,Int},Tdict}
    Wigner6j::ThreadSafeDict{NTuple{6,UInt},Tdict}

    # these arrays replace the global caches in the previous WignerSymbols.jl.
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
    
    # BigInts for numerators and denominators of s ⋅ √r
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

"""
Convenience function for generating a StaticPrimeFactorization of length `nfact`.
"""
factorbuffer(nfact) = StaticPrimeFactorization(zeros(UInt32, nfact), one(Int8),1)

"""
Bounds on the prime-counting function.

Rosser, J. Barkley; Schoenfeld, Lowell (1962). "Approximate formulas for some functions of 
prime numbers". Illinois J. Math. 6: 64–94.
"""
πbound(n) = (n < 2) ? 0 : Int(ceil(1.25506n / log(n)))


"""
This constructs the BoundedWignerCache.

Bounds on the maximum factorial and the number of series terms are based on Table 3 of 
Johansson and Forssén 2016.

Set nj=3 for 3j, nj=6 for 6j.
"""
function BoundedWignerCache(
        dict3j::ThreadSafeDict{Tuple{UInt,UInt,UInt,Int,Int},Tdict}, 
        dict6j::ThreadSafeDict{NTuple{6,UInt},Tdict},
        max_j::Integer, nj=6) where {Tdict}

    # compute the bounds
    multiplier = Int(ceil(2 + round(nj / 3)))
    maxfactorial = multiplier * max_j + 1
    maxt = max_j + 1
    max_prime_factors = πbound(maxfactorial)

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

        Ref{BigInt}(big(0)),
        Ref{BigInt}(big(0)),
        Ref{BigInt}(big(0)),
        
        Ref{BigInt}(big(0)),
        Ref{BigInt}(big(0)),
        Ref{BigInt}(big(0)),
        Ref{BigInt}(big(0)),

        factorbuffer(max_prime_factors),  # s2n
        factorbuffer(max_prime_factors),  # snum 
        factorbuffer(max_prime_factors),  # rnum
        factorbuffer(max_prime_factors),  # sden
        factorbuffer(max_prime_factors)   # rden
    )
    return cache
end

"""
Expand a StaticPrimeFactorization to a new size.
"""
function grow_factorization!(spf::StaticPrimeFactorization{T}, new_size::Integer) where T
    old_size = length(spf.powers)
    resize!(spf.powers, new_size)
    @inbounds for i in old_size:new_size
        spf.powers[i] = zero(T)
    end
end

"""
Change the maximum jᵢ of a BoundedWignerCache. This requires reallocating everything, it's 
a real mess. In `wigner3j` and `wigner6j` we grow the BoundedWignerCache by a factor of 2 
as necessary, i.e. amortized constant time.
"""
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


"""
Computes the prime factors of an integer, and stores that in the cache.
"""
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
        fact_dest::StaticPrimeFactorization{T}, n::Integer) where T
    precompute_primefactorial!(cache, n)
    row = cache.factorialtable[n+1]
    for i in 1:length(row)
        fact_dest.powers[i] += row[i]
    end
    fact_dest.last_nonzero_index = max(length(row), fact_dest.last_nonzero_index)
    return
end

"""
Deposit the factorization of an integer n to a destination fact_dest.
"""
function primefactorial!(cache::BoundedWignerCache, 
                         fact_dest::StaticPrimeFactorization{T}, n::Integer) where {T, N}
    one!(fact_dest, 1)  # reset factor destination
    mul_primefactorial!(cache, fact_dest, n)
end

"""
Convenience function for multiplying together a bunch of primefactorials in a buffer.
"""
function primefactorial!(cache::BoundedWignerCache, 
                         fact_dest::AbstractPrimeFactorization{T}, 
                         ns::NTuple{N, Integer}) where {T, N}
    one!(fact_dest, 1)  # reset factor destination to +1
    for n in ns
        mul_primefactorial!(cache, fact_dest, n)
    end
    return
end

"""
Sets a prime factorization to 1 by setting all the nonzero powers to zero. One can specify
a sign for ±1.
"""
function one!(fact_dest::StaticPrimeFactorization{T}, sign::Integer) where T
    for i in 1:fact_dest.last_nonzero_index
        fact_dest.powers[i] = zero(T)
    end
    # fill!(fact_dest.powers, zero(T))
    fact_dest.sign = sign < 0 ? -one(Int8) : one(Int8)
    fact_dest.last_nonzero_index = 1
    return
end


"""
# Style Note
Our operations which manipulate StaticPrimeFactorization don't ever decrease the last 
nonzero index. It's maybe possible that something cancelled, but tests show that such an
optimization has a negligible effect on performance. We do increase it as necessary.
"""

"""
Divide common factors of StaticPrimeFactorization numbers.
"""
function divgcd!(cache::BoundedWignerCache, 
        a::StaticPrimeFactorization, b::StaticPrimeFactorization)

    af, bf = a.powers, b.powers
    @inbounds for k = 1:min(a.last_nonzero_index, b.last_nonzero_index)
        gk = min(af[k], bf[k])
        af[k] -= gk
        bf[k] -= gk
    end
    return
end

"""
Find the minimum of each power between two factorizations, and put it into the first 
argument `a`.
"""
function _vmin!(a::StaticPrimeFactorization{U}, 
    b::StaticPrimeFactorization{U}) where {U<:Unsigned}
    af = a.powers
    bf = b.powers
    @inbounds for k = 1:a.last_nonzero_index
        af[k] = min(af[k], bf[k])
    end
    return 
end

"""
Subtract the powers of two StaticPrimeFactorization numbers.

The resulting number has the max of the two arguments' last nonzero indices.
"""
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

"""
Add the powers of two StaticPrimeFactorization numbers.
"""
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

"""
Recover a BigInteger from the factorization. This sets the input `A` to the result.
"""
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

"""
Recover a BigInteger from the factorization. This creates a new BigInt.
"""
function _convert(cache::BoundedWignerCache, T::Type{BigInt}, a::StaticPrimeFactorization)
    A = one(BigInt)
    return _convert!(cache, A, a)
end

"""
Auxiliary function to compute sums of a list of PrimeFactorizations as quickly as possible.

This mutates the input vector of StaticPrimeFactorization and returns the sum. This 
unfortunately incurs some allocations by operating recursively.
"""
function sumlist!(cache::BoundedWignerCache, list::Vector{StaticPrimeFactorization{T}}, 
                  ind = 1:length(list)) where {T}
    # first compute gcd to take out common factors
    g = cache.snum
    one!(g, 1)
    first_p = list[ind[1]]
    @inbounds for i in 1:first_p.last_nonzero_index
        g.powers[i] = first_p.powers[i]
    end
    g.last_nonzero_index = first_p.last_nonzero_index
    @inbounds for k in ind
        _vmin!(g, list[k])
    end
    @inbounds for k in ind
        _vsub!(list[k], g)
    end
    gint = _convert(cache, BigInt, g)  # TODO: make this non-allocating
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


"""
Squared triangle coefficient.

This puts the results in the arguments s1n and s1d.
"""
function Δ²!(cache::BoundedWignerCache, s1n::T, s1d::T, 
             j₁, j₂, j₃) where T <: StaticPrimeFactorization
    # also checks the triangle conditions by converting to unsigned integer:
    primefactorial!(cache, s1n, (
        convert(UInt, + j₁ + j₂ - j₃),
        convert(UInt, + j₁ - j₂ + j₃),
        convert(UInt, - j₁ + j₂ + j₃)))

    primefactorial!(cache, s1d, convert(UInt, j₁ + j₂ + j₃ + 1) )
end

"""
Compute the sum appearing in the 3j symbol.

This makes heavy use of preallocated variables in the BoundedWignerCache.
"""
function compute3jseries(cache::BoundedWignerCache, β₁, β₂, β₃, α₁, α₂)
    krange = max(α₁, α₂, zero(α₁)):min(β₁, β₂, β₃)
    numk = length(krange)
    T_int = eltype(eltype(cache.factorialtable))
    T = PrimeFactorization{T_int}

    # get the preallocated variables from the cache
    nums = cache.nums
    dens = cache.dens

    for (i, k) in enumerate(krange)
        num = nums[i]
        den = dens[i]
        one!(num, iseven(k) ? 1 : -1)
        primefactorial!(cache, den, (k, k-α₁, k-α₂, β₁-k, β₂-k, β₃-k))
        # divgcd!(cache, num, den)
    end

     # write gcd to buffer
    commondenominator!(cache, nums, dens, numk) 
    totalden = _convert!(cache, cache.seriesden[], cache.denbuf)
    totalnum = sumlist!(cache, nums, 1:numk)
    return totalnum//totalden
end

"""
Remove one factor from each odd exponent to move into square root.
"""
function splitsquare!(s::T, r::T, a::T) where {T <: StaticPrimeFactorization}
    one!(r, a.sign)
    for i in 1:a.last_nonzero_index
        r.powers[i] = convert(UInt8, isodd(a.powers[i]))
    end

    one!(s, 1)
    for i in 1:a.last_nonzero_index
        s.powers[i] = (a.powers[i])>>1
    end

    r.last_nonzero_index = a.last_nonzero_index
    s.last_nonzero_index = a.last_nonzero_index
end


"""
Compute the Wigner-3j symbol. 

The cache parametric type `Tdict` sets the type that the cache stores results. The type 
`T` is what the final result will be converted to. This enables the default behavior 
of storing results as pairs of r and s where r, s are Rational{BigInt} but offering 
arbitrary output.
"""
function wigner3j(cache::BoundedWignerCache{Tdict}, T::Type{<:Real},
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
        grow!(2 * new_max_j, cache)
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
        # get buffered variables 
        s1n, s1d, s2n = cache.numbuf, cache.denbuf, cache.s2n
        snum, rnum, sden, rden = cache.snum, cache.rnum, cache.sden, cache.rden

        # mutating versions of the functions in the main WignerSymbols.jl file
        Δ²!(cache, s1n, s1d, j₁, j₂, j₃)
        splitsquare!(sden, rden, s1d)
        primefactorial!(cache, s2n, (β₂, β₁ - α₁, β₁ - α₂, β₃, β₃ - α₁, β₂ - α₂))
        s2n_mul_s1n = _vadd!(s1n, s2n)  # multiply s2n and s1n
        splitsquare!(snum, rnum, s2n_mul_s1n)  # split square and store in snum and rnum

        s = _convert!(cache, cache.snumint[], snum) // _convert!(cache, cache.sdenint[], sden)
        r = _convert!(cache, cache.rnumint[], rnum) // _convert!(cache, cache.rdenint[], rden)

        series = compute3jseries(cache, β₁, β₂, β₃, α₁, α₂)

        # multiply s by series
        Base.GMP.MPZ.mul!(s.num, series.num)
        Base.GMP.MPZ.mul!(s.den, series.den)

        # store the r and s pair into the dictionary. 
        rs = convert(Tdict, RSPair(r, s))
        cache.Wigner3j[(β₁, β₂, β₃, α₁, α₂)] = store_result(rs)  # deepcopy if necessary
    end
    return convert(T, sgn)*convert(T, rs)
end



# 6j stuff ---

function wigner6j(cache::BoundedWignerCache{Tdict}, T::Type{<:Real}, 
                  j₁, j₂, j₃, j₄, j₅, j₆) where Tdict <: Number
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

    new_max_j = Int(ceil(max(abs(j₁), abs(j₂), abs(j₃))))
    if(cache.max_j[] < new_max_j)
        grow!(2 * new_max_j, cache)
    end

    # dictionary lookup or compute
    if haskey(cache.Wigner6j, (β₁, β₂, β₃, α₁, α₂, α₃))
        rs = cache.Wigner6j[(β₁, β₂, β₃, α₁, α₂, α₃)]
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

        rs = convert(Tdict, RSPair(r, s))
        cache.Wigner6j[(β₁, β₂, β₃, α₁, α₂, α₃)] = store_result(rs)
    end
    return convert(T, rs.s) * convert(T, signedroot(rs.r))
end


# compute the sum appearing in the 6j symbol
function compute6jseries(cache::BoundedWignerCache, β₁, β₂, β₃, α₁, α₂, α₃, α₄)
                         krange = max(α₁, α₂, α₃, α₄):min(β₁, β₂, β₃)
                         T = PrimeFactorization{eltype(eltype(cache.factorialtable))}
    numk = length(krange)
    nums = cache.nums
    dens = cache.dens

    for (i, k) in enumerate(krange)
        num = nums[i]
        den = dens[i]
        primefactorial!(cache, num, k+1)
        num.sign = Int8(iseven(k) ? 1 : -1)
        primefactorial!(cache, den, (k-α₁, k-α₂, k-α₃, k-α₄, β₁-k, β₂-k, β₃-k))
        divgcd!(cache, num, den)
    end

    commondenominator!(cache, nums, dens, numk) 
    totalden = _convert!(cache, cache.seriesden[], cache.denbuf)
    totalnum = sumlist!(cache, nums, 1:numk)

    return totalnum//totalden
end
