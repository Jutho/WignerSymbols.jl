using WignerSymbols

function test_serial(n::Int)
    new_cache = WignerCache()  # create a new cache
    sum = 0.0
    for i in 1000:n
        sum = sum + wigner3j(new_cache, BigFloat, 5000, i, 5000, 0, 0)
    end
    sum
end

using TimerOutputs
WignerSymbols.__init__()
reset_timer!(WignerSymbols.to)
println(test_serial(1500))
print(WignerSymbols.to)

##  NEW =======
GC.gc()
using WignerSymbols

function test_bounded(n::Int)
    new_cache = BoundedWignerCache(5000)  # create a new cache
    sum = 0.0
    for i in 1000:n
        sum = sum + wigner3j(new_cache, BigFloat, 5000, i, 5000, 0, 0)
    end
    sum
end

using TimerOutputs
WignerSymbols.__init__()
reset_timer!(WignerSymbols.to)
println(test_bounded(1500))
print(WignerSymbols.to)

##
