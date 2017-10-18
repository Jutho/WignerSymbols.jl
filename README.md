# WignerSymbols

[![Build Status](https://travis-ci.org/Jutho/WignerSymbols.jl.svg?branch=master)](https://travis-ci.org/Jutho/WignerSymbols.jl)
[![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](LICENSE.md)
[![Coverage Status](https://coveralls.io/repos/Jutho/WignerSymbols.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/Jutho/WignerSymbols.jl?branch=master)
[![codecov.io](http://codecov.io/github/Jutho/WignerSymbols.jl/coverage.svg?branch=master)](http://codecov.io/github/Jutho/WignerSymbols.jl?branch=master)

Compute Wigner's 3j and 6j symbols, and related quantities such as Clebsch-Gordan coefficients and Racah's symbols.

## Requirements
~~This requires a recent master edition of Julia (i.e. v0.7.0-DEV), because it depends on some changes in `Base.GMP`. In particular, it uses the mutating functions for reducing allocation overhead while working with `BigInts` (namely [JuliaLang/julia#21654](https://github.com/JuliaLang/julia/pull/21654)). It also depends on `Primes.jl` for generating prime numbers.~~

Now works from Julia v0.6.0 onwards (by including/backporting the relevant code from Julia Base in this package).


## Installation
~~Until it is registered, install via `Pkg.clone("https://github.com/Jutho/WignerSymbols.jl.git")`.~~

Now that it is register, install via `Pkg.add("WignerSymbols")`.

## Available functions
While the following function signatures are probably self-explanatory, you can query help for them in the Julia REPL to get further details.
*   `wigner3j(T::Type{<:AbstractFloat} = Float64, j₁, j₂, j₃, m₁, m₂, m₃ = -m₂-m₁) -> ::T`
*   `wigner6j(T::Type{<:AbstractFloat} = Float64, j₁, j₂, j₃, j₄, j₅, j₆) -> ::T`
*   `clebschgordan(T::Type{<:AbstractFloat} = Float64, j₁, j₂, j₃, m₁, m₂, m₃ = m₁+m₂) -> ::T`
*   `racahV(T::Type{<:AbstractFloat} = Float64, j₁, j₂, j₃, m₁, m₂, m₃ = -m₁-m₂) -> ::T`
*   `racahW(T::Type{<:AbstractFloat} = Float64, j₁, j₂, J, j₃, J₁₂, J₂₃) -> ::T`
*   `δ(j₁, j₂, j₃) -> ::Bool`
*   `Δ(T::Type{<:AbstractFloat} = Float64, j₁, j₂, j₃) -> ::T`

## Implementation
Largely based on reading the paper (but not the code):

[1] [H. T. Johansson and C. Forssén, SIAM Journal on Scientific Compututing 38 (2016) 376-384](https://doi.org/10.1137/15M1021908) ([arXiv:1504.08329](https://arxiv.org/abs/1504.08329))

with some additional modifications to further improve efficiency for large `j` (angular momenta quantum numbers).

In particular, 3j and 6j symbols are computed exactly, in the format `√(r) * s` where `r` and `s` are exactly computed as `Rational{BigInt}`,
using an intermediate representation based on prime number factorization. As a consequence thereof, all of the above functions can be called
requesting `BigFloat` precision for the result. There is currently no convenient syntax for obtaining `r` and `s` directly (see TODO).

Most intermediate calculations (prime factorizations of numbers and their factorials, conversion between prime powers and `BigInt`s) are cached to improve the efficiency, but this can result in large use of memory when querying Wigner symbols for large values of `j`.

Also uses ideas from

[2] [J. Rasch and A. C. H. Yu, SIAM Journal on Scientific Compututing 25 (2003), 1416–1428](https://doi.org/10.1137/S1064827503422932)

for caching the computed 3j and 6j symbols.

## Benchmark

## Todo
*   Wigner 9-j symbols, as explained in [1] and based on

    [3] [L. Wei, New formula for 9-j symbols and their direct calculation, Computers in Physics, 12 (1998), 632–634.](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.481.5946&rep=rep1&type=pdf)

*   Convenient syntax to get the exact results in the `√(r) * s` format, but in such a way that it can be used by
    the Julia type system and can be converted afterwards.
