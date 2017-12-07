using Primes.isprime
import Base.divgcd

const primetable = [2,3,5]
const factortable = [UInt8[], UInt8[1], UInt8[0,1], UInt8[2], UInt8[0,0,1]]
const factorialtable = [UInt32[], UInt32[], UInt32[1], UInt32[1,1], UInt32[3,1], UInt32[3,1,1]]
const bigprimetable = [[big(2)], [big(3)], [big(5)]]
const bigone = Ref{BigInt}(big(1))

# Make a prime iterator
struct PrimeIterator
end
primes() = PrimeIterator()

Base.iteratorsize(::Type{PrimeIterator}) = Base.IsInfinite()
Base.iteratoreltype(::Type{PrimeIterator}) = Base.HasEltype()
Base.eltype(::PrimeIterator) = Int

# Get the `n`th prime; store all primes up to the `n`th if not yet available
function prime(n::Int)
    p = last(primetable)
    while length(primetable) < n
        p = p + 2
        while !isprime(p)
            p += 2
        end
        push!(primetable, p)
        push!(bigprimetable, [big(p)])
    end
    @inbounds return primetable[n]
end

Base.start(::PrimeIterator) = 1
Base.next(::PrimeIterator, n) = prime(n), n+1
Base.done(::PrimeIterator, n) = false

# get primes and their powers as `BigInt`, also cache all results
function bigprime(n::Integer, e::Integer=1)
    e == 0 && return bigone[]
    p = prime(n) # triggers computation of prime(n) if necessary
    @inbounds l = length(bigprimetable[n])
    @inbounds while l < e
        # compute next prime power as approximate square of existing results
        k = (l+1)>>1
        push!(bigprimetable[n], bigprimetable[n][k]*bigprimetable[n][l+1-k])
        l += 1
    end
    @inbounds return bigprimetable[n][e]
end

# A custom `Integer` subtype to store an integer as its prime factorization
struct PrimeFactorization{T<:Unsigned} <: Integer
    powers::Vector{T}
    sign::Int8
end
PrimeFactorization(powers::Vector{T}) where {T<:Unsigned} = PrimeFactorization{T}(powers, one(Int8))

# define our own factor function, returning an instance of PrimeFactorization
function primefactor(n::Integer)
    iszero(n) && return PrimeFactorization(UInt8[], zero(Int8))
    sn = n < 0 ? -one(Int8) : one(Int8)
    n = abs(n)
    m = length(factortable)
    while m < abs(n)
        m += 1
        powers = UInt8[] # should be sufficient for all integers up to 2^255
        a = m
        for p in primes()
            f = 0
            anext, r = divrem(a, p)
            while r == 0
                f += 1
                a = anext
                anext, r = divrem(a, p)
            end
            push!(powers, f)
            a == 1 && break
        end
        push!(factortable, powers)
    end
    @inbounds return PrimeFactorization(copy(factortable[n]), sn)
end

function primefactorial(n::Integer)
    n < 0 && throw(DomainError(n))
    m = length(factorialtable)-1
    @inbounds while m < n
        prevfactorial = factorialtable[m+1]
        m += 1
        f = primefactor(m).powers
        powers = copy(prevfactorial)
        if length(f) > length(powers) # can at most be 1 larger
            push!(powers, 0)
        end
        for k = 1:length(f)
            powers[k] += f[k]
        end
        push!(factorialtable, powers)
    end
    @inbounds return PrimeFactorization(copy(factorialtable[n+1]))
end

# Methods for PrimeFactorization:
Base.copy(a::PrimeFactorization) = PrimeFactorization(copy(a.powers), a.sign)

Base.one(::Type{PrimeFactorization{T}}) where {T<:Unsigned} = PrimeFactorization(Vector{T}())
Base.zero(::Type{PrimeFactorization{T}}) where {T<:Unsigned} = PrimeFactorization(Vector{T}(), zero(Int8))

Base.promote_rule(P::Type{<:PrimeFactorization},::Type{<:Integer}) = P
Base.promote_rule(P::Type{<:PrimeFactorization},::Type{BigInt}) = BigInt

Base.convert(P::Type{<:PrimeFactorization}, n::Integer) = convert(P, primefactor(n))
function Base.convert(::Type{BigInt}, a::PrimeFactorization)
    A = big(a.sign)
    for (n, e) in enumerate(a.powers)
        if !iszero(e)
            MPZ.mul!(A, bigprime(n, e))
        end
    end
    return A
end
Base.convert(::Type{PrimeFactorization{T}}, a::PrimeFactorization{T}) where {T<:Unsigned} = a
Base.convert(::Type{PrimeFactorization{T1}}, a::PrimeFactorization{T2}) where {T1<:Unsigned, T2<:Unsigned} = PrimeFactorization(map(T1, a.powers), a.sign)

Base.:(==)(a::PrimeFactorization, b::PrimeFactorization) = a.powers == b.powers && a.sign == b.sign
function Base.:<(a::PrimeFactorization, b::PrimeFactorization)
    if a.sign != b.sign
        return a.sign < b.sign
    elseif a.sign < 0
        return <(-b, -a)
    else
        ag, bg = divgcd(a, b)
        if length(ag.powers) <= length(bg.powers) && all(k->ag.powers[k]<bg.powers[k], 1:length(ag.powers))
            return true
        else
            return convert(BigInt, ag) < convert(BigInt, bg)
        end
    end
end

# Methods for PrimeFactorization: only fast multiplication, and lcm and gcd.
# Addition and subtraction will require conversion to BigInt
Base.sign(a::PrimeFactorization) = a.sign

Base.:-(a::PrimeFactorization) = PrimeFactorization(a.powers, -a.sign)
function Base.:*(a::PrimeFactorization{T}, b::PrimeFactorization{T}) where {T}
    if a.sign == 0
        return a
    elseif b.sign ==0
        return b
    else
        return PrimeFactorization(_vadd!(copy(a.powers), b.powers), a.sign*b.sign)
    end
end
function Base.gcd(a::PrimeFactorization{T}, b::PrimeFactorization{T}) where {T}
    if a.sign == 0
        return b
    elseif b.sign ==0
        return a
    else
        return PrimeFactorization(_vmin!(copy(a.powers), b.powers))
    end
end
function Base.lcm(a::PrimeFactorization{T}, b::PrimeFactorization{T}) where {T}
    if a.sign == 0
        return a
    elseif b.sign ==0
        return b
    else
        return PrimeFactorization(_vmax!(copy(a.powers), b.powers))
    end
end
Base.divgcd(a::PrimeFactorization{T}, b::PrimeFactorization{T}) where {T} = divgcd!(copy(a), copy(b))
function divgcd!(a::PrimeFactorization{T}, b::PrimeFactorization{T}) where {T}
    af, bf = a.powers, b.powers
    for k = 1:min(length(af), length(bf))
        gk = min(af[k], bf[k])
        af[k] -= gk
        bf[k] -= gk
    end
    while length(af) > 0 && iszero(last(af))
        pop!(af)
    end
    while length(bf) > 0 && iszero(last(bf))
        pop!(bf)
    end
    return a, b
end

# split `a::PrimeFactorization` into a square `s` and a remainder `r`, such that
# `a = s^2 * r` and the powers in the prime factorization of `r` are zero or one
function splitsquare(a::PrimeFactorization)
    r = PrimeFactorization(map(p->convert(UInt8, isodd(p)), a.powers), a.sign)
    while length(r.powers) > 0 && iszero(last(r.powers))
        pop!(r.powers)
    end
    s = PrimeFactorization(map(p->(p>>1), a.powers))
    while length(s.powers) > 0 && iszero(last(s.powers))
        pop!(s.powers)
    end
    return s, r
end

# given a list of numerators and denominators, compute the common denominator and
# the rescaled numerator after putting all fractions at the same common denominator
function commondenominator!(nums::Vector{P}, dens::Vector{P}) where {P<:PrimeFactorization}
    isempty(nums) && return one(P)
    # accumulate lcm of denominator
    den = PrimeFactorization(copy(dens[1].powers))
    for i = 2:length(dens)
        _vmax!(den.powers, dens[i].powers)
    end
    # rescale numerators
    for i = 1:length(nums)
        _vsub!(_vadd!(nums[i].powers, den.powers), dens[i].powers)
    end
    return den
end

# auxiliary function to compute sums of a list of PrimeFactorizations as quickly as possible
function sumlist!(list::Vector{<:PrimeFactorization}, ind = 1:length(list))
    # first compute gcd to take out common factors
    g = PrimeFactorization(copy(list[ind[1]].powers))
    for k in ind
        _vmin!(g.powers, list[k].powers)
    end
    for k in ind
        _vsub!(list[k].powers, g.powers)
    end
    L = length(ind)
    if L > 32
        l = L >> 1
        s = sumlist!(list, first(ind).+(0:l-1)) + sumlist!(list, first(ind).+(l:L-1))
    else
        # do sum
        s = big(0)
        for k in ind
            MPZ.add!(s, convert(BigInt, list[k]))
        end
    end
    return MPZ.mul!(s, convert(BigInt, g))
end

# Mutating vector methods that also grow and shrink as required
function _vmin!(af::Vector{U}, bf::Vector{U}) where {U<:Unsigned}
    while length(af) > length(bf)
        pop!(af)
    end
    @inbounds for k = 1:length(af)
        af[k] = min(af[k], bf[k])
    end
    while length(af) > 0 && iszero(last(af))
        pop!(af)
    end
    return af
end
function _vmax!(af::Vector{U}, bf::Vector{U}) where {U<:Unsigned}
    while length(bf) > length(af)
        push!(af, zero(U))
    end
    @inbounds for k = 1:length(bf)
        af[k] = max(af[k], bf[k])
    end
    return af
end
function _vadd!(af::Vector{U}, bf::Vector{U}) where {U<:Unsigned}
    while length(bf) > length(af)
        push!(af, zero(U))
    end
    @inbounds for k = 1:length(bf)
        af[k] = +(af[k], bf[k])
    end
    return af
end
function _vsub!(af::Vector{U}, bf::Vector{U}) where {U<:Unsigned}
    if length(bf) > length(af)
        throw(OverflowError())
    end
    @inbounds for k = 1:length(bf)
        bf[k] > af[k] && throw(OverflowError())
        af[k] -= bf[k]
    end
    while length(af) > 0 && iszero(last(af))
        pop!(af)
    end
    return af
end
