

struct WignerCache
    Wigner3j::Dict{Tuple{UInt,UInt,UInt,Int,Int},Tuple{Rational{BigInt},Rational{BigInt}}}
    Wigner6j::Dict{NTuple{6,UInt},Tuple{Rational{BigInt},Rational{BigInt}}}

    primetable::Array{Int64,1}
    factortable::Array{Array{UInt8,1},1}
    factorialtable::Array{Array{UInt32,1},1}
    bigprimetable::Array{Array{BigInt,1},1}
    bigone::Base.RefValue{BigInt}
end

function WignerCache()
    cache = WignerCache(
        Dict{Tuple{UInt,UInt,UInt,Int,Int},Tuple{Rational{BigInt},Rational{BigInt}}}(),
        Dict{NTuple{6,UInt},Tuple{Rational{BigInt},Rational{BigInt}}}(),
        [2,3,5],
        [UInt8[], UInt8[1], UInt8[0,1], UInt8[2], UInt8[0,0,1]],
        [UInt32[], UInt32[], UInt32[1], UInt32[1,1], UInt32[3,1], UInt32[3,1,1]],
        [[big(2)], [big(3)], [big(5)]],
        Ref{BigInt}(big(1))
    )
    cache.bigone[] = big(1)
    cache.bigprimetable[1][1] = big(2)
    cache.bigprimetable[2][1] = big(3)
    cache.bigprimetable[3][1] = big(5)
    return cache
end

default_cache = WignerCache()

