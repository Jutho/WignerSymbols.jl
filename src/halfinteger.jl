# HalfInteger
struct HalfInteger <: Real
    num::Int
end
Base.:+(a::HalfInteger, b::HalfInteger) = HalfInteger(a.num+b.num)
Base.:-(a::HalfInteger, b::HalfInteger) = HalfInteger(a.num-b.num)
Base.:-(a::HalfInteger) = HalfInteger(-a.num)
Base.:<=(a::HalfInteger, b::HalfInteger) = a.num <= b.num
Base.:<(a::HalfInteger, b::HalfInteger) = a.num < b.num
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
        throw(InexactError())
    end
end
Base.convert(::Type{HalfInteger}, r::Real) = convert(HalfInteger, convert(Rational, r))
Base.convert(T::Type{<:Integer}, s::HalfInteger) = iseven(s.num) ? convert(T, s.num>>1) : throw(InexactError())
Base.convert(T::Type{<:Rational}, s::HalfInteger) = convert(T, s.num//2)
Base.convert(T::Type{<:Real}, s::HalfInteger) = convert(T, s.num/2)
Base.convert(::Type{HalfInteger}, s::HalfInteger) = s

Base.isinteger(a::HalfInteger) = iseven(a.num)
ishalfinteger(a::HalfInteger) = true
ishalfinteger(a::Integer) = true
ishalfinteger(a::Rational) = a.den == 1 || a.den == 2
ishalfinteger(a::Real) = isinteger(2*a)

converthalfinteger(a::Number) = convert(HalfInteger, a)
