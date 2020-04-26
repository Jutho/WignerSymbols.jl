using WignerSymbols
using ThreadSafeDicts

function test_serial(n::Int)
    new_cache = WignerCache()  # create a new cache
    sum = zero(BigFloat)
    for i in 1:n
        for j in 1:n
            for k in 1:n
                sum += wigner3j(new_cache, BigFloat, i, j, k, 0, 0)
            end
        end
    end
    sum
end

function test_bounded(n::Int)
    WignerSymbols.__init__()
    sum = zero(BigFloat)
    for i in 1:n
        for j in 1:n
            for k in 1:n
                sum += wigner3j(BigFloat, i, j, k, 0, 0)
            end
        end
    end
    sum
end


function test_threaded(n::Int)
    WignerSymbols.__init__()
    summands = zeros(BigFloat, Threads.nthreads())
    Threads.@threads for i in 1:n
        tid = Threads.threadid()
        for j in 1:n
            for k in 1:n
                summands[tid] += wigner3j(BigFloat, i, j, k, 0, 0)
            end
        end
    end
    sum(summands)
end

function test_threaded_specialized(n::Int)
    caches = WignerSymbols.get_thread_caches(BigFloat, n)
    summands = zeros(BigFloat, Threads.nthreads())
    Threads.@threads for i in 1:n
        tid = Threads.threadid()
        for j in 1:n
            for k in 1:n
                summands[tid] += wigner3j(caches[tid], BigFloat, i, j, k, 0, 0)
            end
        end
    end
    sum(summands)
end

maxj = 150
##

@time begin
    println("default serial:  ", test_serial(maxj))
    GC.gc()
end
##
@time begin
    println("bounded serial:  ", test_bounded(maxj))
    GC.gc()
end
##
@time begin
    println("bounded thread:  ", test_threaded(maxj))
    GC.gc()
end
##
@time begin
    println("specialized th:  ", test_threaded_specialized(maxj))
    GC.gc()
end
##
