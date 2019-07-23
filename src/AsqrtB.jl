const IntorRational = Union{Integer,Rational}

"""
    AsqrtB{T<:Rational} <: Number

Number type of the form A * sqrt(B), where A and B are rational numbers
"""
struct AsqrtB{T<:Rational} <: Real
    a::T
    b::T
end

AsqrtB(a::IntorRational, b::IntorRational) = AsqrtB(promote(Rational(a),Rational(b))...)

AsqrtB(p) = AsqrtB(p,one(p))
AsqrtB(x::AsqrtB) = x
AsqrtB{T}(p) where T = AsqrtB{T}(p,one(p))
AsqrtB{T}(x::AsqrtB) where T = AsqrtB{T}(x.a,x.b) 

function Base.show(io::IO, x::AsqrtB)
    a,b = x.a,x.b
    if denominator(a)==1 
        a=Int(a)
        a_print = "$a"
    else
        a_print = "$(a.num)/$(a.den)"
    end
    if denominator(b)==1 
        b=Int(b)
        b_print = "$b"
    else
        b_print = "$(b.num)/$(b.den)"
    end
    print(io, a_print," * sqrt(",b_print, ")")
end

Base.:*(x::IntorRational, y::AsqrtB) = AsqrtB(x*y.a, y.b)
Base.:*(x::AsqrtB, y::IntorRational) = *(y,x)
Base.:*(x::AsqrtB, y::AsqrtB) = AsqrtB(x.a*y.a, x.b*y.b)
Base.:-(x::AsqrtB) = AsqrtB(-x.a, x.b)

(T::Type{<:AbstractFloat})(x::AsqrtB) = T(x.a) * sqrt(T(x.b))# T(x.a*sqrt(x.b))
Base.promote_rule(::Type{AsqrtB{T}}, ::Type{AsqrtB{S}}) where {T<:Rational, S<:Rational} = AsqrtB{promote_type(T,S)}
Base.promote_rule(::Type{AsqrtB{T}}, ::Type{S}) where {T<:Rational, S<:Number} = promote_type(T,float(S))

   
