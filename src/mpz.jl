# taken from Julia 0.7-Dev base library

module MPZ
# wrapping of libgmp functions
# - "output parameters" are labeled x, y, z, and are returned when appropriate
# - constant input parameters are labeled a, b, c
# - a method modifying its input has a "!" appendend to its name, according to Julia's conventions
# - some convenient methods are added (in addition to the pure MPZ ones), e.g. `add(a, b) = add!(BigInt(), a, b)`
#   and `add!(x, a) = add!(x, x, a)`.
using Base.GMP: BigInt, Limb

const mpz_t = Ref{BigInt}
const bitcnt_t = Culong

gmpz(op::Symbol) = (Symbol(:__gmpz_, op), :libgmp)

init!(x::BigInt) = (ccall((:__gmpz_init, :libgmp), Void, (mpz_t,), x); x)
init2!(x::BigInt, a) = (ccall((:__gmpz_init2, :libgmp), Void, (mpz_t, bitcnt_t), x, a); x)

realloc2!(x, a) = (ccall((:__gmpz_realloc2, :libgmp), Void, (mpz_t, bitcnt_t), x, a); x)
realloc2(a) = realloc2!(BigInt(), a)

sizeinbase(a::BigInt, b) = Int(ccall((:__gmpz_sizeinbase, :libgmp), Csize_t, (mpz_t, Cint), a, b))

for op in (:add, :sub, :mul, :fdiv_q, :tdiv_q, :fdiv_r, :tdiv_r, :gcd, :lcm, :and, :ior, :xor)
    op! = Symbol(op, :!)
    @eval begin
        $op!(x::BigInt, a::BigInt, b::BigInt) = (ccall($(gmpz(op)), Void, (mpz_t, mpz_t, mpz_t), x, a, b); x)
        $op(a::BigInt, b::BigInt) = $op!(BigInt(), a, b)
        $op!(x::BigInt, b::BigInt) = $op!(x, x, b)
    end
end

invert!(x::BigInt, a::BigInt, b::BigInt) =
    ccall((:__gmpz_invert, :libgmp), Cint, (mpz_t, mpz_t, mpz_t), x, a, b)
invert(a::BigInt, b::BigInt) = invert!(BigInt(), a, b)
invert!(x::BigInt, b::BigInt) = invert!(x, x, b)

for op in (:add_ui, :sub_ui, :mul_ui, :mul_2exp, :fdiv_q_2exp, :pow_ui, :bin_ui)
    op! = Symbol(op, :!)
    @eval begin
        $op!(x::BigInt, a::BigInt, b) = (ccall($(gmpz(op)), Void, (mpz_t, mpz_t, Culong), x, a, b); x)
        $op(a::BigInt, b) = $op!(BigInt(), a, b)
        $op!(x::BigInt, b) = $op!(x, x, b)
    end
end

ui_sub!(x::BigInt, a, b::BigInt) = (ccall((:__gmpz_ui_sub, :libgmp), Void, (mpz_t, Culong, mpz_t), x, a, b); x)
ui_sub(a, b::BigInt) = ui_sub!(BigInt(), a, b)

for op in (:scan1, :scan0)
    @eval $op(a::BigInt, b) = Int(ccall($(gmpz(op)), Culong, (mpz_t, Culong), a, b))
end

mul_si!(x::BigInt, a::BigInt, b) = (ccall((:__gmpz_mul_si, :libgmp), Void, (mpz_t, mpz_t, Clong), x, a, b); x)
mul_si(a::BigInt, b) = mul_si!(BigInt(), a, b)
mul_si!(x::BigInt, b) = mul_si!(x, x, b)

for op in (:neg, :com, :sqrt, :set)
    op! = Symbol(op, :!)
    @eval begin
        $op!(x::BigInt, a::BigInt) = (ccall($(gmpz(op)), Void, (mpz_t, mpz_t), x, a); x)
        $op(a::BigInt) = $op!(BigInt(), a)
    end
    op == :set && continue # MPZ.set!(x) would make no sense
    @eval $op!(x::BigInt) = $op!(x, x)
end

for (op, T) in ((:fac_ui, Culong), (:set_ui, Culong), (:set_si, Clong), (:set_d, Cdouble))
    op! = Symbol(op, :!)
    @eval begin
        $op!(x::BigInt, a) = (ccall($(gmpz(op)), Void, (mpz_t, $T), x, a); x)
        $op(a) = $op!(BigInt(), a)
    end
end

popcount(a::BigInt) = Int(ccall((:__gmpz_popcount, :libgmp), Culong, (mpz_t,), a))

mpn_popcount(d::Ptr{Limb}, s::Integer) = Int(ccall((:__gmpn_popcount, :libgmp), Culong, (Ptr{Limb}, Csize_t), d, s))
mpn_popcount(a::BigInt) = mpn_popcount(a.d, abs(a.size))

function tdiv_qr!(x::BigInt, y::BigInt, a::BigInt, b::BigInt)
    ccall((:__gmpz_tdiv_qr, :libgmp), Void, (mpz_t, mpz_t, mpz_t, mpz_t), x, y, a, b)
    x, y
end
tdiv_qr(a::BigInt, b::BigInt) = tdiv_qr!(BigInt(), BigInt(), a, b)

powm!(x::BigInt, a::BigInt, b::BigInt, c::BigInt) =
    (ccall((:__gmpz_powm, :libgmp), Void, (mpz_t, mpz_t, mpz_t, mpz_t), x, a, b, c); x)
powm(a::BigInt, b::BigInt, c::BigInt) = powm!(BigInt(), a, b, c)
powm!(x::BigInt, b::BigInt, c::BigInt) = powm!(x, x, b, c)

function gcdext!(x::BigInt, y::BigInt, z::BigInt, a::BigInt, b::BigInt)
    ccall((:__gmpz_gcdext, :libgmp), Void, (mpz_t, mpz_t, mpz_t, mpz_t, mpz_t), x, y, z, a, b)
    x, y, z
end
gcdext(a::BigInt, b::BigInt) = gcdext!(BigInt(), BigInt(), BigInt(), a, b)

cmp(a::BigInt, b::BigInt) = Int(ccall((:__gmpz_cmp, :libgmp), Cint, (mpz_t, mpz_t), a, b))
cmp_si(a::BigInt, b) = Int(ccall((:__gmpz_cmp_si, :libgmp), Cint, (mpz_t, Clong), a, b))
cmp_ui(a::BigInt, b) = Int(ccall((:__gmpz_cmp_ui, :libgmp), Cint, (mpz_t, Culong), a, b))
cmp_d(a::BigInt, b) = Int(ccall((:__gmpz_cmp_d, :libgmp), Cint, (mpz_t, Cdouble), a, b))

mpn_cmp(a::Ptr{Limb}, b::Ptr{Limb}, c) = ccall((:__gmpn_cmp, :libgmp), Cint, (Ptr{Limb}, Ptr{Limb}, Clong), a, b, c)
mpn_cmp(a::BigInt, b::BigInt, c) = mpn_cmp(a.d, b.d, c)

get_str!(x, a, b::BigInt) = (ccall((:__gmpz_get_str,:libgmp), Ptr{Cchar}, (Ptr{Cchar}, Cint, mpz_t), x, a, b); x)
set_str!(x::BigInt, a, b) = Int(ccall((:__gmpz_set_str, :libgmp), Cint, (mpz_t, Ptr{UInt8}, Cint), x, a, b))
get_d(a::BigInt) = ccall((:__gmpz_get_d, :libgmp), Cdouble, (mpz_t,), a)

limbs_write!(x::BigInt, a) = ccall((:__gmpz_limbs_write, :libgmp), Ptr{Limb}, (mpz_t, Clong), x, a)
limbs_finish!(x::BigInt, a) = ccall((:__gmpz_limbs_finish, :libgmp), Void, (mpz_t, Clong), x, a)
import!(x::BigInt, a, b, c, d, e, f) = ccall((:__gmpz_import, :libgmp), Void,
    (mpz_t, Csize_t, Cint, Csize_t, Cint, Csize_t, Ptr{Void}), x, a, b, c, d, e, f)

end # module MPZ
