using Test
using WignerSymbols

## wigner3j summation tests
function test_serial_3j(n::Int)
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

function test_threaded_3j(n::Int)
    WignerSymbols.__init__()  # reset global caches
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

@testset "wigner3j threaded sum" begin
    for max_j in [7, 30, 100]
        tol = 10*max_j*eps(BigFloat)
        @test test_serial_3j(max_j) ≈ test_threaded_3j(max_j) atol=tol
    end
end

## wigner6j summation tests
function test_serial_6j(n::Int)
    new_cache = WignerCache()  # create a new cache
    s = zero(BigFloat)
    for i in 1:n
        s += wigner6j(new_cache, BigFloat, i, i, i, i, i, i)
    end
    s
end

function test_threaded_6j(n::Int)
    summands = zeros(BigFloat, Threads.nthreads())
    Threads.@threads for i in 1:n
        summands[Threads.threadid()] += wigner6j(BigFloat, i, i, i, i, i, i)
    end
    sum(summands)
end

@testset "wigner6j threaded sum" begin
    for max_j in [7, 51, 100, 300, 701]
        tol = 10*max_j*eps(BigFloat)
        @test test_serial_6j(max_j) ≈ test_threaded_6j(max_j) atol=tol
    end
end
