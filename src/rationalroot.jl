#module Tmp
using Primes: factor
export RationalRoot

const IntorRational = Union{Integer,Rational}

"""
    RationalRoot{T<:Rational} <: Real

Signed square root of a rational number
"""
struct RationalRoot{T<:Rational} <: Real
    #Representation of the number sign(r)*sqrt(|r|) where r is rational
    r::T
end

RationalRoot(n::Integer) = RationalRoot(Rational(n))
RationalRoot(x::RationalRoot) = x
RationalRoot{T}(x::RationalRoot) where {T<:Rational} = RationalRoot{T}(x.r) 

Base.:-(x::RationalRoot) = RationalRoot(-x.r)
Base.:*(x::IntorRational, y::RationalRoot) = RationalRoot(y.r*x*abs(x))
Base.:*(x::RationalRoot, y::IntorRational) = *(y,x)
Base.:*(x::RationalRoot, y::RationalRoot) = RationalRoot(x.r*y.r)
Base.:/(x::RationalRoot, y::RationalRoot) = RationalRoot(x.r/y.r)
Base.:/(x::IntorRational, y::RationalRoot) = x*inv(y)
Base.:/(x::RationalRoot, y::IntorRational) = inv(y/x)
Base.:+(x::RationalRoot, y::RationalRoot) = float(x) + float(y)
Base.:-(x::RationalRoot, y::RationalRoot) = float(x) - float(y)

"""
    simplify(x::RationalRoot)

return rational numbers (a,b) such that x = a * sqrt(b)
"""
function simplify(x::RationalRoot)
    r = x.r
    sgn = sign(r)
    r *= sgn    # ensure r is positive from here
    num, den = r.num, r.den
    num_f = factor(num)
    den_f = factor(den)

    num_a=num_b=1
    for (p,i) in num_f
        num_a *= p^(div(i,2))
        num_b *= p^(i%2)
    end
    den_a=den_b=1
    for (p,i) in den_f
        den_a *= p^(div(i,2))
        den_b *= p^(i%2)
    end
    return (sgn*num_a//den_a, num_b//den_b) 
end

function Base.show(io::IO, x::RationalRoot)
    a,b = simplify(x)
    if a==0 || b==0
        return print(io, "0")
    end
    if denominator(a)==1 
        a=Int(a)
        a_print = "$a"
    else
        a_print = "$(a.num)/$(a.den)"
    end
    if b==1
        return print(io, a_print)
    end
    if denominator(b)==1 
        b=Int(b)
        b_print = "$b"
    else
        b_print = "$(b.num)/$(b.den)"
    end
    if a==1
        return print(io,  "sqrt(",b_print,")")
    elseif a==-1
        return print(io, "-sqrt(",b_print,")")
    end
    print(io, a_print," * sqrt(",b_print, ")")
end

# Conversion to float
(T::Type{<:AbstractFloat})(x::RationalRoot) = sign(x.r) * sqrt(T(abs(x.r)))

Base.promote_rule(::Type{RationalRoot{T}}, ::Type{RationalRoot{S}}) where {T<:Rational, S<:Rational} = RationalRoot{promote_type(T,S)}
Base.promote_rule(::Type{RationalRoot{T}}, ::Type{S}) where {T<:Rational, S<:Number} = promote_type(T,float(S))

#end
