# HalfInteger
"""
    struct HalfInteger <: Real

Represents half-integer values.
"""
struct HalfInteger <: Real
    numerator::Int # with an implicit denominator of 2

    """
        HalfInteger(numerator::Integer, denominator::Integer)

    Constructs a `HalfInteger` object as a rational number from the given integer numerator
    and denominator values.
    """
    function HalfInteger(num::Integer, den::Integer)
        (den == 2) && return new(num)
        (den == 1) && return new(2*num)
        (den == 0) && throw(ArgumentError("Denominator can not be zero."))
        # If non-trivial, we'll see if we can reduce it down to a half-integer
        numerator, r = divrem(2*num, den)
        if r == 0
            return new(numerator)
        else
            throw(ArgumentError("$num // $den is not a half-integer value."))
        end
    end
end
HalfInteger(x::Real) = convert(HalfInteger, x)

Base.promote_rule(::Type{HalfInteger}, ::Type{<:Integer}) = HalfInteger
Base.promote_rule(::Type{HalfInteger}, T::Type{<:Rational}) = T
Base.promote_rule(::Type{HalfInteger}, T::Type{<:Real}) = T

Base.convert(::Type{HalfInteger}, n::Integer) = HalfInteger(2*n, 2)
function Base.convert(::Type{HalfInteger}, r::Rational)
    if r.den == 1
        return HalfInteger(2*r.num, 2)
    elseif r.den == 2
        return HalfInteger(r.num, 2)
    else
        throw(InexactError(:HalfInteger, HalfInteger, r))
    end
end
function Base.convert(::Type{HalfInteger}, r::Real)
    num = 2*r
    if isinteger(num)
        return HalfInteger(convert(Int, num), 2)
    else
        throw(InexactError(:HalfInteger, HalfInteger, r))
    end
end
Base.convert(T::Type{<:Integer}, s::HalfInteger) = iseven(s.numerator) ? convert(T, s.numerator>>1) : throw(InexactError(Symbol(T), T, s))
Base.convert(T::Type{<:Rational}, s::HalfInteger) = convert(T, s.numerator//2)
Base.convert(T::Type{<:Real}, s::HalfInteger) = convert(T, s.numerator/2)
Base.convert(::Type{HalfInteger}, s::HalfInteger) = s

# Arithmetic

Base.:+(a::HalfInteger, b::HalfInteger) = HalfInteger(a.numerator+b.numerator, 2)
Base.:-(a::HalfInteger, b::HalfInteger) = HalfInteger(a.numerator-b.numerator, 2)
Base.:-(a::HalfInteger) = HalfInteger(-a.numerator, 2)
Base.:*(a::Integer, b::HalfInteger) = HalfInteger(a * b.numerator, 2)
Base.:*(a::HalfInteger, b::Integer) = b * a
Base.:<=(a::HalfInteger, b::HalfInteger) = a.numerator <= b.numerator
Base.:<(a::HalfInteger, b::HalfInteger) = a.numerator < b.numerator
Base.one(::Type{HalfInteger}) = HalfInteger(2, 2)
Base.zero(::Type{HalfInteger}) = HalfInteger(0, 2)

# Hashing

function Base.hash(a::HalfInteger, h::UInt)
    iseven(a.numerator) && return hash(a.numerator>>1, h)
    num, den = a.numerator, 2
    den = 1
    pow = -1
    if abs(num) < 9007199254740992
        return hash(ldexp(Float64(num),pow), h)
    end
    h = Base.hash_integer(den, h)
    h = Base.hash_integer(pow, h)
    h = Base.hash_integer(num, h)
    return h
end

# Parsing and printing

"""
    parse(HalfInteger, s)

Parses the string `s` into the corresponding `HalfInteger`-value. String can either be a
number or a fraction of the form `<x>/2`.
"""
function Base.parse(::Type{HalfInteger}, s::AbstractString)
    if in('/', s)
        num, den = split(s, '/'; limit=2)
        parse(Int, den) == 2 ||
            throw(ArgumentError("Denominator not 2 in HalfInteger string '$s'."))
        HalfInteger(parse(Int, num), 2)
    elseif !isempty(strip(s))
        HalfInteger(parse(Int, s))
    else
        throw(ArgumentError("input string is empty or only contains whitespace"))
    end
end

Base.show(io::IO, x::HalfInteger) =
    print(io, iseven(x.numerator) ? "$(div(x.numerator, 2))" : "$(x.numerator)/2")

# Other methods

Base.isinteger(a::HalfInteger) = iseven(a.numerator)
ishalfinteger(a::HalfInteger) = true
ishalfinteger(a::Integer) = true
ishalfinteger(a::Rational) = a.den == 1 || a.den == 2
ishalfinteger(a::Real) = isinteger(2*a)

converthalfinteger(a::Number) = convert(HalfInteger, a)

Base.numerator(a::HalfInteger) = iseven(a.numerator) ? div(a.numerator, 2) : a.numerator
Base.denominator(a::HalfInteger) = iseven(a.numerator) ? 1 : 2

# Range of HalfIntegers

struct HalfIntegerRange <: AbstractVector{HalfInteger}
    start :: HalfInteger
    stop :: HalfInteger

    function HalfIntegerRange(start::HalfInteger, stop::HalfInteger)
        (start <= stop) ||
            throw(ArgumentError("Second argument must be greater or equal to the first."))
        isinteger(stop - start) ||
            throw(ArgumentError("Two arguments must have integer difference."))
        return new(start, stop)
    end
end
Base.iterate(it::HalfIntegerRange) = (it.start, it.start + 1)
Base.iterate(it::HalfIntegerRange, s) = (s <= it.stop) ? (s, s+1) : nothing
Base.length(it::HalfIntegerRange) = convert(Int, it.stop - it.start) + 1
Base.size(it::HalfIntegerRange) = (length(it),)
function Base.getindex(it::HalfIntegerRange, i::Integer)
    1 <= i <= length(it) || throw(BoundsError(it, i))
    it.start + i - 1
end

Base.:(:)(i::HalfInteger, j::HalfInteger) = HalfIntegerRange(i, j)
