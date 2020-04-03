using WignerSymbols

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
    new_cache = BoundedWignerCache(n)  # create a new cache
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

function setup_caches(jmax)
    caches = BoundedWignerCache[]
    resize!(empty!(caches), Threads.nthreads())
    for i in 1:Threads.nthreads()
        cache = BoundedWignerCache(jmax)
        caches[i] = cache
    end
    return caches
end

function test_threaded(n::Int, caches)
    summands = zeros(BigFloat, Threads.nthreads())
    Threads.@threads for i in 1:n
        for j in 1:n
            for k in 1:n
                threadid = Threads.threadid()
                summands[threadid] += wigner3j(caches[threadid], BigFloat, i, j, k, 0, 0)
            end
        end
    end
    sum(summands)
end

maxj = 120
@time println("standard       ", test_serial(maxj))
GC.gc()
@time println("bounded cache: ", test_bounded(maxj))
GC.gc()

caches = setup_caches(maxj)
@time println("threaded:      ", test_threaded(maxj, caches))
GC.gc()

caches = setup_caches(maxj)
GC.enable(false)
@time println("thread no GC:  ", test_threaded(maxj, caches))
GC.enable(true)
GC.gc()
##
