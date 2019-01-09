using Test
using WignerSymbols: HalfInteger, ishalfinteger, HalfIntegerRange

@testset "HalfInteger" begin
    @testset "HalfInteger type" begin
        # HalfInteger constructors
        @test HalfInteger(1, 2).numerator == 1
        @test HalfInteger(1, 1).numerator == 2
        @test HalfInteger(0, 1).numerator == 0
        @test HalfInteger(0, 2).numerator == 0
        @test HalfInteger(0, 5).numerator == 0
        @test HalfInteger(10, 5).numerator == 4
        @test HalfInteger(21, 14).numerator == 3
        @test HalfInteger(-3, 2).numerator == -3
        @test HalfInteger(3, -2).numerator == -3
        @test HalfInteger(-3, -2).numerator == 3
        @test_throws ArgumentError HalfInteger(1, 0)
        @test_throws ArgumentError HalfInteger(1, 3)
        @test_throws ArgumentError HalfInteger(1, -3)
        @test_throws ArgumentError HalfInteger(-5, 3)
        @test_throws ArgumentError HalfInteger(-1000, -999)

        # convert methods
        @test convert(HalfInteger, 2) == HalfInteger(2, 1)
        @test convert(HalfInteger, 1//2) == HalfInteger(1, 2)
        @test convert(HalfInteger, 1.5) == HalfInteger(3, 2)
        @test_throws InexactError convert(HalfInteger, 1//3)
        @test_throws InexactError convert(HalfInteger, 0.6)
        @test convert(HalfInteger, 2) == 2
        @test convert(HalfInteger, 1//2) == 1//2
        @test convert(HalfInteger, 1.5) == 1.5
        @test_throws InexactError convert(Integer, HalfInteger(1, 2))

        # single-argument constructor
        @test HalfInteger(0) == HalfInteger(0, 2)
        @test HalfInteger(1) == HalfInteger(1, 1)
        @test HalfInteger(2) == HalfInteger(2, 1)
        @test HalfInteger(-30) == HalfInteger(-60, 2)
        @test HalfInteger(0//2) == HalfInteger(0, 1)
        @test HalfInteger(1//2) == HalfInteger(1, 2)
        @test HalfInteger(-5//2) == HalfInteger(-5, 2)
    end

    a = HalfInteger(2)
    b = HalfInteger(3, 2)

    @testset "HalfInteger arithmetic" begin
        @test a + b == 2 + 3//2
        @test a - b == 2 - 3//2
        @test zero(a) == 0
        @test one(a) == 1
        @test a > b
        @test b < a
        @test b <= a
        @test a >= b
        @test a == a
        @test a != b
        @test 2 * HalfInteger(0) == HalfInteger(0)
        @test 2 * HalfInteger(1, 2) == HalfInteger(1)
        @test HalfInteger(1) * 2 == HalfInteger(2)
        @test 2 * a == HalfInteger(4)
        @test (-1) * b == HalfInteger(-3//2)

        @test floor(HalfInteger(0)) === HalfInteger(0)
        @test floor(HalfInteger(1, 2)) === HalfInteger(0)
        @test floor(HalfInteger(-1, 2)) === HalfInteger(-1)
        @test floor(HalfInteger(-1)) === HalfInteger(-1)
        @test floor(HalfInteger(1)) === HalfInteger(1)
        @test floor(HalfInteger(5, 2)) === HalfInteger(2)
        @test floor(HalfInteger(-5, 2)) === HalfInteger(-3)
        @test floor(HalfInteger(-5)) === HalfInteger(-5)
        @test floor(HalfInteger(5)) === HalfInteger(5)

        @test floor(Int, HalfInteger(0)) === 0
        @test floor(Int, HalfInteger(1, 2)) === 0
        @test floor(Int32, HalfInteger(-5, 2)) === Int32(-3)
        @test floor(Int32, HalfInteger(5)) === Int32(5)
    end

    @testset "Parsing and printing" begin
        @test string(HalfInteger(0)) == "0"
        @test string(HalfInteger(1)) == "1"
        @test string(HalfInteger(-1)) == "-1"
        @test string(HalfInteger(1, 2)) == "1/2"
        @test string(HalfInteger(-3, 2)) == "-3/2"

        # hi"" string macro / parse(::HalfInteger)
        @test parse(HalfInteger, "0") == HalfInteger(0)
        @test parse(HalfInteger, "1") == HalfInteger(1)
        @test parse(HalfInteger, "210938") == HalfInteger(210938)
        @test parse(HalfInteger, "-15") == HalfInteger(-15)
        @test parse(HalfInteger, "1/2") == HalfInteger(1//2)
        @test parse(HalfInteger, "-3/2") == HalfInteger(-3//2)
        @test_throws ArgumentError parse(HalfInteger, "")
        @test_throws ArgumentError parse(HalfInteger, "-50/100")
        @test_throws ArgumentError parse(HalfInteger, "1/3")
    end

    @testset "HalfInteger hashing" begin
        @test hash(a) == hash(2)
        @test hash(b) == hash(1.5)
    end

    @testset "Other HalfInteger methods" begin
        @test isinteger(HalfInteger(0))
        @test isinteger(HalfInteger(1))
        @test !isinteger(HalfInteger(1, 2))

        @test ishalfinteger(1)
        @test ishalfinteger(1.0)
        @test ishalfinteger(-0.5)
        @test ishalfinteger(HalfInteger(0))
        @test ishalfinteger(HalfInteger(1, 2))
        @test ishalfinteger(1//1)
        @test ishalfinteger(1//2)
        @test !ishalfinteger(0.3)
        @test !ishalfinteger(-5//7)

        @test numerator(HalfInteger(0)) == 0
        @test numerator(HalfInteger(1, 2)) == 1
        @test numerator(HalfInteger(1)) == 1
        @test numerator(HalfInteger(-3, 2)) == -3

        @test denominator(HalfInteger(0)) == 1
        @test denominator(HalfInteger(1, 2)) == 2
        @test denominator(HalfInteger(1)) == 1
        @test denominator(HalfInteger(-3, 2)) == 2
    end

    @testset "HalfIntegerRange" begin
        hi(x) = HalfInteger(x)

        @test length(HalfIntegerRange(hi(0), hi(0))) == 1
        @test length(HalfIntegerRange(hi(0), hi(2))) == 3
        let hirange = HalfIntegerRange(hi(-1//2), hi(1//2))
            @test length(hirange) == 2
            @test size(hirange) == (2,)
            @test collect(hirange) == [hi(-1//2), hi(1//2)]
        end
        let hirange = HalfIntegerRange(hi(0), hi(1//2))
            @test length(hirange) == 1
            @test size(hirange) == (1,)
            @test collect(hirange) == [hi(0)]
        end
        let hirange = HalfIntegerRange(hi(1//2), hi(3))
            @test length(hirange) == 3
            @test size(hirange) == (3,)
            @test collect(hirange) == [hi(1//2), hi(3//2), hi(5//2)]
        end

        @test hi(5):hi(7) == HalfIntegerRange(hi(5), hi(7))
        @test hi(-1//2):hi(1//2) == HalfIntegerRange(hi(-1//2), hi(1//2))

        @test collect(hi(0) : hi(2)) == [hi(0), hi(1), hi(2)]
        @test collect(hi(-3//2) : hi(1//2)) == [hi(-3//2), hi(-1//2), hi(1//2)]

        let hirange = hi(-3//2):hi(0)
            @test length(hirange) == 2
            @test size(hirange) == (2,)
            @test collect(hirange) == [hi(-3//2), hi(-1//2)]
        end

        @test hi(1//2) ∈ hi(-1//2) : hi(1//2)
        @test 1 ∈ hi(0) : hi(2)
        @test 1//2 ∈ hi(-1//2) : hi(7//2)
        @test !(hi(1//2) ∈ hi(0) : hi(1))
        @test !(1//2 ∈ hi(-1) : hi(7))

        r = hi(-3//2) : hi(3//2)
        @test r[1] == hi(-3//2)
        @test r[2] == hi(-1//2)
        @test r[3] == hi(1//2)
        @test r[4] == hi(3//2)
        @test_throws BoundsError r[0]
        @test_throws BoundsError r[5]
    end
end
