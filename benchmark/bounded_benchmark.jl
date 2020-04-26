using WignerSymbols
using ThreadSafeDicts

function test_serial(n::Int)
    new_cache = WignerCache()  # create a new cache
    sum = 0.0
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
    sum = 0.0
    for i in 1:n
        for j in 1:n
            for k in 1:n
                sum += wigner3j(Float64, i, j, k, 0, 0)
            end
        end
    end
    sum
end


function test_threaded(n::Int)
    WignerSymbols.__init__()
    summands = zeros(Float64, Threads.nthreads())
    Threads.@threads for i in 1:n
        tid = Threads.threadid()
        for j in 1:n
            for k in 1:n
                summands[tid] += wigner3j(Float64, i, j, k, 0, 0)
            end
        end
    end
    sum(summands)
end

function test_threaded_specialized(n::Int)
    caches = WignerSymbols.get_thread_caches(Float64, n)
    summands = zeros(Float64, Threads.nthreads())
    Threads.@threads for i in 1:n
        tid = Threads.threadid()
        for j in 1:n
            for k in 1:n
                summands[tid] += wigner3j(caches[tid], Float64, i, j, k, 0, 0)
            end
        end
    end
    sum(summands)
end

maxj = 150
##
@time println("default serial:  ", test_serial(maxj))
GC.gc()
##
@time println("bounded serial:  ", test_bounded(maxj))
GC.gc()
##
@time println("bounded thread:  ", test_threaded(maxj))
GC.gc()
##
@time println("bounded thread specialized:  ", test_threaded_specialized(maxj))
GC.gc()
##

# caches, mydict = setup_caches(maxj)
# GC.enable(false)
# @time println("thread no GC:    ", test_threaded(maxj, caches, mydict))
# GC.enable(true)
# GC.gc()

##