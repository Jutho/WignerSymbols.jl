#module Tmp
export RationalRoot, rationalsqrt

const IntorRational{T} = Union{T,Rational{T}} where {T<:Integer}
const _inner = Val{:inner}()
"""
    RationalRoot{T} <: Real
Signed square root of a rational number
"""
struct RationalRoot{T<:Integer} <: AbstractIrrational
    #Representation of the number s*sqrt(r)
    s::Rational{T} # signed rational
    r::T # positive integer

    function RationalRoot{T}(::Val{:inner}, s::Rational{T}, r::T) where {T<:Integer}
        r > 0 || throw(DomainError(r, "factor in square root should be positive"))
        r2 = primefactor(r)
        for i = 1:length(r2.powers)
            while r2.powers[i] > 1
                p = prime(i)
                s *= convert(T, p)
                r = div(r, p*p)
                r2.powers[i] -= 2
            end
        end
        return new{T}(s, r)
    end
end

RationalRoot(::Val{:inner}, s::Rational{T}, r::T) where {T<:Integer} =
    RationalRoot{T}(_inner, s, r)

RationalRoot(n::Integer) = RationalRoot(n//one(n))
RationalRoot(s::Rational{T}) where {T} = RationalRoot{T}(_inner, s, one(T))
RationalRoot(x::RationalRoot) = x

RationalRoot{T}(n::Integer) where {T<:Integer }= RationalRoot{T}(convert(T, n)//one(T))
RationalRoot{T}(s::Rational) where {T} =
    RationalRoot{T}(Val{:inner}(), convert(Rational{T}, s), one(T))
RationalRoot{T}(x::RationalRoot) where {T} =
    RationalRoot{T}(Val{:inner}(), convert(Rational{T}, x.s), convert(T, x.r))

rationalsqrt(x::T) where {T<:Integer} = RationalRoot{T}(_inner, one(T)//one(T), x)
rationalsqrt(x::Rational{T}) where {T<:Integer} =
    RationalRoot{T}(_inner, one(T)//x.den, x.num*x.den)

Base.:(==)(a::RationalRoot, b::RationalRoot) = a.s == b.s && a.r == b.r
Base.:(==)(a::RationalRoot, b::Real) = a.r == one(a.r) && a.s == b
Base.:(==)(a::Real, b::RationalRoot) = b == a

Base.one(::Type{RationalRoot{T}}) where T = RationalRoot{T}(one(T))
Base.zero(::Type{RationalRoot{T}}) where T = RationalRoot{T}(zero(T))
Base.isone(x::RationalRoot) = x.r == one(x.r) && isone(x.s)
Base.iszero(x::RationalRoot) = x.r == one(x.r) && iszero(x.s)
Base.isinteger(x::RationalRoot) = x.r == one(x.r) && isinteger(x.s)

function Base.hash(a::RationalRoot, h::UInt)
    h = hash(a.s, h)
    if a.r != one(a.r)
        h = hash(a.r, h)
    end
    return h
end

Base.convert(::Type{T}, x::RationalRoot) where {T<:AbstractFloat} =
    convert(T, x.s.num)/convert(T, x.s.den) * sqrt(convert(T, x.r))
Base.convert(::Type{T}, x::RationalRoot) where {T<:IntorRational} =
    x.r == one(x.r) ? convert(T, x.s) : throw(InexactError(nameof(T), T, x))

Base.AbstractFloat(x::RationalRoot) = convert(AbstractFloat, x)
Base.Float32(x::RationalRoot) = convert(Float32, x)
Base.Float64(x::RationalRoot) = convert(Float64, x)
Base.BigFloat(x::RationalRoot) = convert(BigFloat, x)

function Base.show(io::IO, x::RationalRoot)
    show(io, x.s)
    print(io, " ⋅ √")
    show(io, x.r)
end

Base.:-(x::RationalRoot) = RationalRoot(Val{:inner}(), -x.s, x.r)
Base.:*(x::IntorRational{T}, y::RationalRoot) where {T<:Integer} =
    RationalRoot(_inner, x*y.s, one(T)*y.r)
Base.:*(x::RationalRoot, y::IntorRational) = *(y,x)
Base.:*(x::RationalRoot, y::RationalRoot) = RationalRoot(_inner, x.s*y.s, x.r*y.r)

Base.inv(x::RationalRoot) = RationalRoot(_inner, inv(x.s*x.r), x.r)

Base.:/(x::RationalRoot, y::RationalRoot) = RationalRoot(_inner, x.s/(y.s*y.r) , x.r*y.r)
Base.:/(x::IntorRational, y::RationalRoot) = x*inv(y)
Base.:/(x::RationalRoot, y::IntorRational) = *(one(y)//y, x)

Base.promote_rule(::Type{RationalRoot{T}}, ::Type{RationalRoot{S}}) where
    {T<:Integer, S<:Integer} = RationalRoot{promote_type(T,S)}
