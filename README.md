# BracketAlgebras.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://Saschobolt.github.io/BracketAlgebras.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://Saschobolt.github.io/BracketAlgebras.jl/dev/)
[![Build Status](https://github.com/Saschobolt/BracketAlgebras.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/Saschobolt/BracketAlgebras.jl/actions/workflows/CI.yml?query=branch%3Amaster)

BracketAlgebras.jl is a package to construct the bracket algebra - the ring of $d\times d$-minors of a matrix in $\mathbb{C}^{n\times d}$. Bracket algebras implement the Ring Interface of the package [AbstractAlgebra.jl](https://github.com/Nemocas/AbstractAlgebra.jl). 

## Overview
To use the full functionality of BracketAlgebras.jl, load both AbstractAlgebra and BracketAlgebras
```julia
julia> using AbstractAlgebra

julia> using BracketAlgebras
```

Bracket algebras can be constructed by providing the number of rows and columns (-1) of the matrix, it models.
```julia
julia> B = BracketAlgebra(6,2)
Bracket algebra over P^2 on 6 points with point ordering 1 < 2 < 3 < 4 < 5 < 6 and coefficient ring Integers.
```

Bracket algebra elements can be added, subtracted, multiplied and raised to a power.
```julia
julia> b = (B([2,3,4])^2 * B([1,3,5]) + B([4,5,6]))^2
[2, 3, 4]^4[1, 3, 5]^2 + 2[4, 5, 6][2, 3, 4]^2[1, 3, 5] + [4, 5, 6]^2
```
The normal form of a bracket algebra expression can be computed by applying the straightening algorithm.
```julia
julia> straighten(b)
[2, 3, 5]^2[2, 3, 4]^2[1, 3, 4]^2 - 2[3, 4, 5][2, 3, 5][2, 3, 4]^2[1, 3, 4][1, 2, 3] + [3, 4, 5]^2[2, 3, 4]^2[1, 2, 3]^2 + 2[4, 5, 6][2, 3, 5][2, 3, 4][1, 3, 4] - 2[4, 5, 6][3, 4, 5][2, 3, 4][1, 2, 3] + [4, 5, 6]^2
```

The full documentation is available under https://saschobolt.github.io/BracketAlgebras.jl/
