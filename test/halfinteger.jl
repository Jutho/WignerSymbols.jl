using Test
using WignerSymbols: HalfInteger

@testset "HalfInteger" begin
    @testset "HalfInteger type" begin
        # HalfInteger constructors
        @test HalfInteger(1, 2).twofold == 1
        @test HalfInteger(1, 1).twofold == 2
        @test HalfInteger(0, 1).twofold == 0
        @test HalfInteger(0, 2).twofold == 0
        @test HalfInteger(0, 5).twofold == 0
        @test HalfInteger(10, 5).twofold == 4
        @test HalfInteger(21, 14).twofold == 3
        @test HalfInteger(-3, 2).twofold == -3
        @test HalfInteger(3, -2).twofold == -3
        @test HalfInteger(-3, -2).twofold == 3
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
end
