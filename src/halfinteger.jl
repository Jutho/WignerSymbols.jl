# HalfInteger
struct HalfInteger <: Real
    twofold::Int
end
Base.:+(a::HalfInteger, b::HalfInteger) = HalfInteger(a.twofold+b.twofold)
Base.:-(a::HalfInteger, b::HalfInteger) = HalfInteger(a.twofold-b.twofold)
Base.:-(a::HalfInteger) = HalfInteger(-a.twofold)
Base.:<=(a::HalfInteger, b::HalfInteger) = a.twofold <= b.twofold
Base.:<(a::HalfInteger, b::HalfInteger) = a.twofold < b.twofold
Base.one(::Type{HalfInteger}) = HalfInteger(2)
Base.zero(::Type{HalfInteger}) = HalfInteger(0)

Base.promote_rule(::Type{HalfInteger}, ::Type{<:Integer}) = HalfInteger
Base.promote_rule(::Type{HalfInteger}, T::Type{<:Rational}) = T
Base.promote_rule(::Type{HalfInteger}, T::Type{<:Real}) = T

Base.convert(::Type{HalfInteger}, n::Integer) = HalfInteger(2*n)
function Base.convert(::Type{HalfInteger}, r::Rational)
    if r.den == 1
        return HalfInteger(2*r.num)
    elseif r.den == 2
        return HalfInteger(r.num)
    else
        throw(InexactError(:HalfInteger, HalfInteger, r))
    end
end
function Base.convert(::Type{HalfInteger}, r::Real)
    num = 2*r
    if isinteger(num)
        return HalfInteger(convert(Int, num))
    else
        throw(InexactError(:HalfInteger, HalfInteger, r))
    end
end
Base.convert(T::Type{<:Integer}, s::HalfInteger) = iseven(s.twofold) ? convert(T, s.twofold>>1) : throw(InexactError(Symbol(T), T, s))
Base.convert(T::Type{<:Rational}, s::HalfInteger) = convert(T, s.twofold//2)
Base.convert(T::Type{<:Real}, s::HalfInteger) = convert(T, s.twofold/2)
Base.convert(::Type{HalfInteger}, s::HalfInteger) = s

function Base.hash(a::HalfInteger, h::UInt)
    iseven(a.twofold) && return hash(a.twofold>>1, h)
    num, den = a.twofold, 2
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

Base.isinteger(a::HalfInteger) = iseven(a.twofold)
ishalfinteger(a::HalfInteger) = true
ishalfinteger(a::Integer) = true
ishalfinteger(a::Rational) = a.den == 1 || a.den == 2
ishalfinteger(a::Real) = isinteger(2*a)

converthalfinteger(a::Number) = convert(HalfInteger, a)
